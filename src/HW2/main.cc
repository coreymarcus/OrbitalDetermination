#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Util.h" // utility functions
#include <iomanip>      // std::setprecision

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//namespaces, don't ever use name space, simply call them aliases
// namespace mat = matlab::engine;

int main() {

	//form a propagator object
	VehicleState::Propagator propobj;

	//set physics constants
	propobj.mu_ = 398600.4; // km^3/sec^2
	propobj.J2_ = 0.00108248;
	propobj.Rearth_ = 6378.145; //km
	propobj.earthrotationspeed_ = 7.29211585530066 * pow(10,-5);
	propobj.C_D_ = 2.0;
	propobj.A_ = 3.6; // m^2
	propobj.m_ = 1350; //kg
	propobj.rho_0_ = 4*pow(10,-13); //kg/m^3
	propobj.r0_ = 7298.145; //km
	propobj.H_ = 200.0; //km

	//parameters
	propobj.useJ2_ = true;
	propobj.usedrag_ = false;

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

	//set tolerance options
	propobj.abstol_ = 1*pow(10,-16);
	propobj.reltol_ = 3*pow(10,-14);
	propobj.dt_var_ = 0.1;

	//duplicate object for using drag
	VehicleState::Propagator propobj_drag = propobj;
	propobj_drag.usedrag_ = true;

	//propagate the orbit for two revolutions
	double dt = 20; //seconds for propagation
	const int N = 86400/20; // propagate for one day
	Eigen::MatrixXd OEhist(7,N); //orbital elements
	Eigen::MatrixXd thist(1,N); //time
	Eigen::MatrixXd Ehist(2,N); //energy
	Eigen::MatrixXd hhist(3,N); //angular momentum
	Eigen::MatrixXd OEhist_drag(7,N); //orbital elements
	Eigen::MatrixXd Ehist_drag(2,N); //energy
	for (int ii = 0; ii < N; ++ii){

		//propagate
		propobj.Propagate(dt);
		propobj_drag.Propagate(dt);

		//update orbital elements
		propobj.State2OE();
		propobj_drag.State2OE();

		//store variables
		OEhist(0,ii) = propobj.a_;
		OEhist(1,ii) = propobj.e_;
		OEhist(2,ii) = propobj.i_;
		OEhist(3,ii) = propobj.ascend_;
		OEhist(4,ii) = propobj.periap_;
		OEhist(5,ii) = propobj.T_p_;
		OEhist(6,ii) = propobj.P_;
		thist(0,ii) = propobj.t_;
		Ehist(0,ii) = propobj.GetKEsp();
		Ehist(1,ii) = propobj.GetPEsp();
		hhist.block(0,ii,3,1) = propobj.GetAngMomVector();

		OEhist_drag(0,ii) = propobj_drag.a_;
		OEhist_drag(1,ii) = propobj_drag.e_;
		OEhist_drag(2,ii) = propobj_drag.i_;
		OEhist_drag(3,ii) = propobj_drag.ascend_;
		OEhist_drag(4,ii) = propobj_drag.periap_;
		OEhist_drag(5,ii) = propobj_drag.T_p_;
		OEhist_drag(6,ii) = propobj_drag.P_;
		Ehist_drag(0,ii) = propobj_drag.GetKEsp();
		Ehist_drag(1,ii) = propobj_drag.GetPEsp();
	}

	std::cout << std::setprecision(17);
	std::cout << "Propagation with J2: \n";
	std::cout << "final position:\n" << propobj.pos_ << std::endl;
	std::cout << "final velocity:\n" << propobj.vel_ << std::endl;

	std::cout << "\n" << "Propagation with J2 and Drag: \n";
	std::cout << "final position:\n" << propobj_drag.pos_ << std::endl;
	std::cout << "final velocity:\n" << propobj_drag.vel_ << std::endl;

	//write data to csv
	Util::Eigen2csv("../data/OEhist_HW2.csv",OEhist);
	Util::Eigen2csv("../data/thist_HW2.csv",thist);
	Util::Eigen2csv("../data/Ehist_HW2.csv",Ehist);
	Util::Eigen2csv("../data/hkhist_HW2.csv",hhist);
	Util::Eigen2csv("../data/OEhistDrag_HW2.csv",OEhist_drag);
	Util::Eigen2csv("../data/EhistDrag_HW2.csv",Ehist_drag);

	std::cout << "done!" << std::endl;
	
	
} // main