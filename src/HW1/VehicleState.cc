#include "VehicleState.h"
#include <iostream> //outputs
#include <Eigen/Dense> // vectors
#include <math.h> // trig functions

namespace VehicleState {

	//default constructor for Propagator
	Propagator::Propagator(){}

	void Propagator::OE2State(){

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

	} // Propagate


} //namespace VehicleState