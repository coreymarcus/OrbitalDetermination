#include "VehicleState.h"
#include <iostream> //outputs
#include <Eigen/Dense> // vectors
#include <math.h> // trig functions
#include "Util.h" //utility functions
#include <vector> // basic cpp vectors
#include <boost/numeric/odeint.hpp> //integrator

namespace VehicleState {

	//default constructor for Propagator
	Propagator::Propagator(){}

	void Propagator::OE2State(){

		//extract locals
		double mu = this->mu_; // gravity constant
		double a = this->a_; //semi-major axis
		double e = this->e_; //eccentricity
		double i = this->i_; //inclination
		double Ohm = this->ascend_; //longitude of ascending node
		double w = this->periap_; //argument of periapsis
		double nu = this->nu_; //true anomaly

		//convert true anomaly to radian
		double nu_rad = nu*180.0/M_PI;

		//calculate the norm of position
		double nr = a*(1.0 - pow(e,2)) / (1 + e*cos(nu_rad));

		//angular velocity
		double h = sqrt(mu*nr*(1+e*cos(nu_rad)));

		// radius in perifocal reference frame
		Eigen::Vector3d r_pqw(cos(nu_rad), sin(nu_rad), 0.0);
		r_pqw = pow(h,2)/(mu*(1+e*cos(nu_rad)))*r_pqw;

		// velocity in perifocal reference frame
		Eigen::Vector3d v_pqw(-sin(nu_rad), e+cos(nu_rad), 0);
		v_pqw = v_pqw*(mu/h);

		// generate 313 rotation matric from perifocal to geocentric vector rotation
		Eigen::Vector3d angles(M_PI*Ohm/180, M_PI*i/180, M_PI*w/180);
		std::vector<int> order{3,1,3};
		Eigen::Matrix3d R = Util::Angle2RotM(angles, order);

		// rotate vectors
		Eigen::Vector3d pos = R*r_pqw;
		Eigen::Vector3d vel = R*v_pqw;

		//Display
		std::cout << "position: \n" << pos << std::endl;
		std::cout << "velocity: \n" << vel << std::endl;

		//Assign

	} // OE2State

	void Propagator::State2OE(){

		// extract locals
		double mu = this->mu_;
		Eigen::Vector3d pos = this->pos_;
		Eigen::Vector3d vel = this->vel_;

		//normalized pos and vel
		double nr = pos.norm();
		double nv = vel.norm();

		//find the angular momentum vector
		Eigen::Vector3d h = pos.cross(vel);

		//now we can calculate semi-major axis
		double a = -mu/(pow(nv,2) - 2*mu/nr);

		//ecentricity vector
		Eigen::Vector3d ecen = (1/mu)*(vel.cross(h) - (mu/nr)*pos);

		//ecentricity
		double e = ecen.norm();

		//inclination
		double inc = acos(h[2]/h.norm()) * 180.0 / M_PI;

		// now we need the line of nodes
		Eigen::Vector3d x(0.0, 0.0, 1.0);
		Eigen::Vector3d nodeline = x.cross(h);

		//calculate longitude of ascending node
		double Ohm;
		if (nodeline[1] >= 0) {
			Ohm = acos(nodeline[0]/nodeline.norm()) * 180.0 / M_PI;
		} else {
			Ohm = 360.0 - acos(nodeline[0]/nodeline.norm()) * 180.0 / M_PI;
		}

		// normalized eccentricity, line of nodes, and position
		Eigen::Vector3d nhat = nodeline/nodeline.norm();
		Eigen::Vector3d ehat = ecen/ecen.norm();
		Eigen::Vector3d rhat = pos/pos.norm();

		// calculate argument of periapsis
		double w;
		if (ehat[2] >= 0) {
			w = acos(nhat.dot(ehat)) * 180.0 / M_PI; // Q1 or Q4
		} else {
			w = 360.0 - acos(nhat.dot(ehat)) * 180.0 / M_PI; // Q2 or Q3
		}

		//calculate true anomally
		double theta;
		if (ehat.dot(rhat) >= 0) {
			theta = acos(ehat.dot(rhat)) * 180.0 / M_PI; // Q1 or Q4
		} else {
			theta = 360.0 - acos(ehat.dot(rhat)) * 180.0 / M_PI; // Q2 or Q3
		}

		// output our results
		std::cout << "semi-major axis [km]: " << a << std::endl;
		std::cout << "eccentricity: " << e << std::endl;
		std::cout << "inclination [deg]: " << inc << std::endl;
		std::cout << "longitude of ascending node [deg]: " << Ohm << std::endl;
		std::cout << "argument of periapsis [deg]: " << w << std::endl;
		std::cout << "true anomally [deg]: " << theta << std::endl;

		//assign results to class
		this->a_ = a; //semi-major axis
		this->e_ = e; //eccentricity
		this->i_ = inc; //inclination
		this->ascend_ = Ohm; //longitude of ascending node
		this->periap_ = w; //argument of periapsis
		this->nu_ = theta; //true anomaly

	} // State2OE

	void Propagator::Propagate(double dt){

		//build state
		state_type xint(6,0.0);
		xint[0] = this->pos_[0];
		xint[1] = this->pos_[1];
		xint[2] = this->pos_[2];
		xint[3] = this->vel_[0];
		xint[4] = this->vel_[1];
		xint[5] = this->vel_[2];

		// std::cout << "before: " << this->pos_ << std::endl;

		//times
		double t1 = this->t_;
		double t2 = t1 + dt;

		//call integrator
		ode::integrate_adaptive( ode::make_controlled( this->abstol_ , this->reltol_ ,
			error_stepper_type() ) , *this , xint , t1 , t2 , this->dt_var_ );

		//update class members
		this->pos_[0] = xint[0];
		this->pos_[1] = xint[1];
		this->pos_[2] = xint[2];
		this->vel_[0] = xint[3];
		this->vel_[1] = xint[4];
		this->vel_[2] = xint[5];
		this->t_ = t2;

		// std::cout << "after: "<< this->pos_ << std::endl;

	} // Propagate

	Eigen::Vector3d Propagator::GetAccelVector(){

		//extract locals
		double mu = this->mu_;
		Eigen::Vector3d pos = this->pos_;

		//radius
		double R = sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2));

		//create accel vector
		Eigen::Vector3d accel = -1.0*mu/pow(R,3)*pos;

		//return
		return accel;


	} // GetAccelVector

	Eigen::Vector3d Propagator::GetAngMomVector(){

		//simple
		return this->pos_.cross(this->vel_);

	} // GetAngMomVector

	double Propagator::GetKEsp(){
		return 0.5*pow(this->vel_.norm(), 2);
	} // GetKEsp

	double Propagator::GetPEsp(){
		return -1.0*this->mu_ / this->pos_.norm();
	} // GetPEsp

	void Propagator::operator()  (const state_type &x , state_type &dxdt , const double /* t */ ){
		
		//extract locals
		double mu = this->mu_;
		state_type pos = {x[0], x[1], x[2]};
		state_type vel = {x[3], x[4], x[5]};

		//calculate acceleration
		double R = sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2));
		double k = -mu/pow(R,3);
		state_type accel = {k*x[0], k*x[1], k*x[2]};

		//change in state
		dxdt[0] = vel[0];
		dxdt[1] = vel[1];
		dxdt[2] = vel[2];
		dxdt[3] = accel[0];
		dxdt[4] = accel[1];
		dxdt[5] = accel[2];

	} //operator()


} //namespace VehicleState