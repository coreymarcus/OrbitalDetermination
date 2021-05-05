#include "VehicleState.h"
#include <iostream> //outputs
#include <Eigen/Dense> // vectors
#include <math.h> // trig functions
#include "Util.h" //utility functions
#include <vector> // basic cpp vectors
#include <boost/numeric/odeint.hpp> //integrator

namespace VehicleState {

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
		// std::cout << "position: \n" << pos << std::endl;
		// std::cout << "velocity: \n" << vel << std::endl;

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
		if (vel.dot(pos) >= 0) {
			theta = acos(ehat.dot(rhat)) * 180.0 / M_PI; // Q1 or Q4
		} else {
			theta = 360.0 - acos(ehat.dot(rhat)) * 180.0 / M_PI; // Q2 or Q3
		}

		// convert theta to radians
		double theta_rad = theta*M_PI/180.0;

		//calculate orbital period
		double P = 2*M_PI*sqrt(pow(a,3)/mu);

		// //eccentric anomaly (radians)
		// double E = acos((e + cos(theta_rad))/(1 + e*cos(theta_rad)));

		// //mean anomally (radians)
		// double M = E - e*sin(E);

		// // time of periapsis passage
		// double T_p = this->t_ - M*P/(2*M_PI);

		//find mean motion
		double meanmot = 1.0/sqrt(a*a*a/mu);

		//eccentric annomaly
		double E = 2.0*atan( tan(theta_rad/2.0) * sqrt( (1.0 - e)/(1.0 + e)));

		//guard domain
		if(E < 0.0) {
			E = E + 2.0*M_PI;
		}

		//mean anomally (radians)
		double M = E - e*sin(E);

		// time of periapsis passage
		double T_p = this->t_ - M/meanmot;

		// output our results
		// std::cout << "semi-major axis [km]: " << a << std::endl;
		// std::cout << "eccentricity: " << e << std::endl;
		// std::cout << "inclination [deg]: " << inc << std::endl;
		// std::cout << "longitude of ascending node [deg]: " << Ohm << std::endl;
		// std::cout << "argument of periapsis [deg]: " << w << std::endl;
		// std::cout << "true anomally [deg]: " << theta << std::endl;
		// std::cout << "time of periapsis passage: " T_p << std::endl;

		//assign results to class
		this->a_ = a; //semi-major axis
		this->e_ = e; //eccentricity
		this->i_ = inc; //inclination
		this->ascend_ = Ohm; //longitude of ascending node
		this->periap_ = w; //argument of periapsis
		this->nu_ = theta; //true anomaly
		this->P_ = P; //orbital period 
		this->T_p_ = T_p; //time from periapse

	} // State2OE

	void Propagator::Propagate(double dt, bool intSTM){

		//set flag
		this->intSTM_ = intSTM;

		//build state
		state_type xint;
		if (intSTM) {
			xint.assign(42,0.0);
			xint[0] = this->pos_[0];
			xint[1] = this->pos_[1];
			xint[2] = this->pos_[2];
			xint[3] = this->vel_[0];
			xint[4] = this->vel_[1];
			xint[5] = this->vel_[2];

			//stm initialization
			xint[6] = 1.0;
			xint[13] = 1.0;
			xint[20] = 1.0;
			xint[27] = 1.0;
			xint[34] = 1.0;
			xint[41] = 1.0;

		} else {
			xint.assign(6,0.0);
			xint[0] = this->pos_[0];
			xint[1] = this->pos_[1];
			xint[2] = this->pos_[2];
			xint[3] = this->vel_[0];
			xint[4] = this->vel_[1];
			xint[5] = this->vel_[2];
		} //fi

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

		//if we propagated the STM, write it out
		if(intSTM){
			Eigen::MatrixXd STM = Eigen::MatrixXd::Zero(6,6);

			for (int ii = 0; ii < 6; ++ii) {
				std::vector<double> col(6,0);
				col[0] = xint[6*ii + 0 + 6];
				col[1] = xint[6*ii + 1 + 6];
				col[2] = xint[6*ii + 2 + 6];
				col[3] = xint[6*ii + 3 + 6];
				col[4] = xint[6*ii + 4 + 6];
				col[5] = xint[6*ii + 5 + 6];

				//convert to eigen
				Eigen::VectorXd eigcol = Util::StdVec2Eigen(col);

				//assign
				STM.block(0,ii,6,1) = eigcol;

			} // for

			//update
			this->STM_ = STM;

		} //fi

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

		//locals
		double mu = this->mu_;
		double Rearth = this->Rearth_;
		double J2 = this->J2_;
		Eigen::Vector3d pos = this->pos_;
		double r = pos.norm();

		return mu / r * (1.0 - J2*pow(Rearth/r,2) * (1.5*pow(pos[2]/r,2) - 0.5) );
	} // GetPEsp

	void Propagator::operator()  (const state_type &x , state_type &dxdt , const double t){
		
		//extract locals
		double mu = this->mu_;
		double Rearth = this->Rearth_;
		double J2 = this->J2_;
		double J3 = this->J3_;
		state_type pos = {x[0], x[1], x[2]};
		state_type vel = {x[3], x[4], x[5]};
		double earthrot = this->earthrotationspeed_;
		double C_D = this->C_D_; // coefficient of drag
		double A = this->A_; // effective area, m^2
		double m = this->m_; // vehicle mass, kg
		double rho_0 = this->rho_0_; // standard air density, kg/m^3
		double r0 = this->r0_; // param for air density calc, km
		double H = this->H_; //param for air density calc, km
		double t_JD_init = this->t_JD_; //UTC julian date at initialization
		double AU2km = this->AU_;
		double mu_sun = this->mu_sun_;
		double mu_moon = this->mu_moon_;

		//initialize some variables we may need
		Eigen::Vector3d r_sun_ECI; //position of the sun in ECI frame [km]
		Eigen::Vector3d r_craft2sun; //ECI position from craft to sun [km]
		Eigen::Vector3d r_moon_ECI; //position of the moon in ECI frame [km]
		Eigen::Vector3d r_craft2moon; //ECI position from craft to moon [km]
		Eigen::Vector3d r_craft(x[0], x[1], x[2]); //eigen position vector for spacecraft [km]

		//intermediate calcs
		double r2 = pow(pos[0],2.0) + pow(pos[1],2.0) + pow(pos[2],2.0);
		double r = sqrt(r2);
		double JD_UTC = t_JD_init + t/(24.0*60.0*60.0); //current UTC Julian Date in days
		double d2r = M_PI/180.0; //degrees to radians

		// std::cout << "t: " << t << " state: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5]; 

		//initalize acceleration
		state_type accel = {0.0, 0.0, 0.0};

		if(this->use20x20_){

			Eigen::Vector3d accelinit = this->gravmodel_->GetGravAccel(r_craft, JD_UTC);

			accel[0] = accel[0] + accelinit[0];
			accel[1] = accel[1] + accelinit[1];
			accel[2] = accel[2] + accelinit[2];

		} else {

			//two body acceleration
			double k = -mu/pow(r,3);
			accel[0] = accel[0] + k*pos[0];
			accel[1] = accel[1] + k*pos[1];
			accel[2] = accel[2] + k*pos[2];

			if (this->useJ2_) { // J-2 accel
				
				double t2 = pow(pos[0],2);
				double t3 = pow(pos[1],2);
				double t4 = pow(pos[2],2);
				double accel_j2_x = J2*pow(Rearth,2)*mu*pos[0]*(t2+t3-t4*4.0)*1.0/pow(t2+t3+t4,7.0/2.0)*(-3.0/2.0);
				double accel_j2_y = J2*pow(Rearth,2)*mu*pos[1]*(t2+t3-t4*4.0)*1.0/pow(t2+t3+t4,7.0/2.0)*(-3.0/2.0);
				double accel_j2_z = J2*pow(Rearth,2)*mu*pos[2]*1.0/pow(t2+t3+t4,7.0/2.0)*(t2*3.0+t3*3.0-t4*2.0)*(-3.0/2.0);

				// add J2 acceleration
				accel[0] = accel[0] + accel_j2_x;
				accel[1] = accel[1] + accel_j2_y;
				accel[2] = accel[2] + accel_j2_z;
			} // fi

		}	

		if (this->usedrag_) { // drag acceleration 
			
			//relative wind vector
			Eigen::Vector3d V_A;
			V_A[0] = vel[0] + earthrot*pos[1];
			V_A[1] = vel[1] - earthrot*pos[0];
			V_A[2] = vel[2];
			double nV_A = V_A.norm(); //magnitude

			//current air density
			double rho_A = rho_0*exp(-(r - r0)/H);

			//add drag acceleration (extra 1000 for unit conversion)
			accel[0] = accel[0] - 1000.0*0.5*C_D*A*rho_A*nV_A*V_A[0]/m;
			accel[1] = accel[1] - 1000.0*0.5*C_D*A*rho_A*nV_A*V_A[1]/m;
			accel[2] = accel[2] - 1000.0*0.5*C_D*A*rho_A*nV_A*V_A[2]/m;
		} //fi

		//first, find the position of the sun in the ECI frame (vallado algo 29)

		//assume JD_UTC = JD_UT1 (less than 1 second deviation)
		double T_UT1 = (JD_UTC - 2451545.0)/36525.0; // julian centuries
		// double T_UT1 = (2453827.5 - 2451545.0)/36525.0; // vallado numbers for checking
		double lambda_M_sun = 280.460 + 36000.771285*T_UT1; // mean sun longitude [degrees]

		//assume T_TBD = T_UT1
		double M_sun = 357.528 + 35999.050957*T_UT1; // mean sun something [degrees]
		double lambda_ecliptic = lambda_M_sun + 1.914666471*sin(M_sun*d2r) + 0.019994643*sin(2*M_sun*d2r);

		//norm of position to sun [AU]
		double r_sun_ECI_norm = 1.000140612 - 0.016708617*cos(M_sun*d2r) - 0.000139589*cos(2*M_sun*d2r);
		double eps = 23.439291 - 0.0130042*T_UT1; //ecentricity of sun orbit?

		//position of sun
		r_sun_ECI[0] = r_sun_ECI_norm*cos(lambda_ecliptic*d2r)*AU2km;
		r_sun_ECI[1] = r_sun_ECI_norm*cos(eps*d2r)*sin(lambda_ecliptic*d2r)*AU2km;
		r_sun_ECI[2] = r_sun_ECI_norm*sin(eps*d2r)*sin(lambda_ecliptic*d2r)*AU2km;

		//relative direction of sun in ECI frame
		r_craft2sun = r_sun_ECI - r_craft;

		bool in_sun = true; //boolean to track if we are in the sun or not

		if (this->useSRP_) { //SRP acceleration

			//check to see if we are in the earth's shadow (vallado algo 34)
			double sundotcraft = r_sun_ECI.dot(r_craft);
			if (sundotcraft < 0)	{ // if positive we are between sun and earth

				//constants for earth shadow
				double alpha_pen = 0.269007205*d2r;
				
				//find the angle between the two [radians]
				double xi = acos(-1.0*sundotcraft/(r_sun_ECI.norm()*r_craft.norm()));

				//satellite horizontal and vertical position
				double sat_horiz = r_craft.norm()*cos(xi);
				double sat_vert = r_craft.norm()*sin(xi);

				double x = Rearth/sin(alpha_pen);

				double pen_vert = tan(alpha_pen)*(x+sat_horiz);
				if(sat_vert <= pen_vert) {
					in_sun = false; //consider only penumbra
				} //fi

			} //fi

			//if we are in the sun, add the SRP force
			if (in_sun)	{
				
				//approximate Area * C_s
				// double SRPcoeff = 0.04*15.0 + 0.59*7.0;
				double SRPcoeff = 1.8*15.0;

				//vallado srp
				double p_srp = 4.57*pow(10.0,-6.0);

				//normalize
				Eigen::Vector3d r_craft2sun_norm = r_craft2sun/r_craft2sun.norm();

				//add SRP acceleration (extra 1000 for unit conversion)
				accel[0] = accel[0] - SRPcoeff*p_srp*r_craft2sun_norm[0]/(1000.0*m);
				accel[1] = accel[1] - SRPcoeff*p_srp*r_craft2sun_norm[1]/(1000.0*m);
				accel[2] = accel[2] - SRPcoeff*p_srp*r_craft2sun_norm[2]/(1000.0*m);
			}


		} //fi

		if (this->useLuniSolar_) { //gravitational effects of the sun and moon

			//first, find the position of the moon in the ECI frame (vallado algo 31)

			//assume JD_UTC = JD_TBD
			// double T_TBD = (2449470.5 - 2451545.0)/36525.0; // vallado numbers
			double T_TBD = (JD_UTC - 2451545.0)/36525.0; // julian centuries
			
			double lambda_ecliptic = 218.32 + 481267.8813*T_TBD + 6.29*sin(d2r*(134.9 + 477198.85*T_TBD))
				- 1.27*sin(d2r*(259.2 - 413335.38*T_TBD)) + 0.66*sin(d2r*(235.7 + 890543.23*T_TBD))
				+ 0.21*sin(d2r*(269.9+954397.70*T_TBD)) - 0.19*sin(d2r*(357.5 + 35999.05*T_TBD))
				- 0.11*sin(186.6 + 966404.05*T_TBD);

			double phi_ecliptic = 5.13*sin(d2r*(93.3 + 483202.03*T_TBD)) + 0.28*sin(d2r*(228.2 + 960400.87*T_TBD))
				- 0.28*sin(d2r*(318.3 + 6003.18*T_TBD)) - 0.17*sin(d2r*(217.6 - 407332.20*T_TBD));

			double D = 0.9508 + 0.0518*cos(d2r*(134.9 + 477198.85*T_TBD)) + 0.0095*cos(d2r*(259.2 - 413335.38*T_TBD))
				+ 0.0078*cos(d2r*(235.7 + 890534.23*T_TBD)) + 0.0028*cos(d2r*(269.9 + 954397.70*T_TBD));

			double eps = 23.439291 - 0.0130042*T_TBD - 1.64*pow(10.0,-7.0)*pow(T_TBD,2.0) + 5.04*pow(10.0,-7.0)*pow(T_TBD,3.0);

			double r_moon_ECI_norm = Rearth/sin(d2r*D);

			r_moon_ECI[0] = r_moon_ECI_norm*cos(d2r*phi_ecliptic)*cos(d2r*lambda_ecliptic);
			r_moon_ECI[1] = r_moon_ECI_norm*(cos(eps*d2r)*cos(d2r*phi_ecliptic)*sin(d2r*lambda_ecliptic) - sin(d2r*eps)*sin(d2r*phi_ecliptic));
			r_moon_ECI[2] = r_moon_ECI_norm*(sin(eps*d2r)*cos(d2r*phi_ecliptic)*sin(d2r*lambda_ecliptic) + cos(d2r*eps)*sin(d2r*phi_ecliptic));

			//relative direction of moon in ECI frame
			r_craft2moon = r_moon_ECI - r_craft;

			//intermediate calcs
			double r2 = pow(pos[0],2.0) + pow(pos[1],2.0) + pow(pos[2],2.0);
			double r = sqrt(r2);
			double r3craft2sun = pow(r_craft2sun.norm(),3.0);
			double r3earth2sun = pow(AU2km*r_sun_ECI_norm,3.0);
			double r3craft2moon = pow(r_craft2moon.norm(),3.0);
			double r3earth2moon = pow(r_moon_ECI_norm,3.0);

			// std::cout << "\n" << accel[0] << ", " << accel[1] << ", " << accel[2] << "\n";

			//two body acceleration to sun
			accel[0] = accel[0] + mu_sun*(r_craft2sun[0]/r3craft2sun - r_sun_ECI[0]/r3earth2sun);
			accel[1] = accel[1] + mu_sun*(r_craft2sun[1]/r3craft2sun - r_sun_ECI[1]/r3earth2sun);
			accel[2] = accel[2] + mu_sun*(r_craft2sun[2]/r3craft2sun - r_sun_ECI[2]/r3earth2sun);

			// std::cout << "\n" << accel[0] << ", " << accel[1] << ", " << accel[2] << "\n";

			//two body acceleration to moon
			accel[0] = accel[0] + mu_moon*(r_craft2moon[0]/r3craft2moon - r_moon_ECI[0]/r3earth2moon);
			accel[1] = accel[1] + mu_moon*(r_craft2moon[1]/r3craft2moon - r_moon_ECI[1]/r3earth2moon);
			accel[2] = accel[2] + mu_moon*(r_craft2moon[2]/r3craft2moon - r_moon_ECI[2]/r3earth2moon);

			// std::cout << "\n" << accel[0] << ", " << accel[1] << ", " << accel[2] << "\n";
			// exit(0);
			
		}

		//integrate the STM if needed
		if (this->intSTM_) {

			//current air density
			double rho_A = rho_0*exp(-(r - r0)/H);
			double drag_coeff = 1000.0*1000.0*0.5*C_D*A*rho_A/m;

			//get the jacobian
			Eigen::Vector3d eigpos(pos.data());
			Eigen::Vector3d eigvel(vel.data());
			Eigen::MatrixXd jac = Util::GetGravJac(eigpos, eigvel, Rearth, earthrot, mu, J2, J3, drag_coeff);

			//extract the STM
			Eigen::MatrixXd STM = Eigen::MatrixXd::Zero(6,6);

			for (int ii = 0; ii < 6; ++ii) {

				std::vector<double> col(6,0);
				col[0] = x[6*ii + 0 + 6];
				col[1] = x[6*ii + 1 + 6];
				col[2] = x[6*ii + 2 + 6];
				col[3] = x[6*ii + 3 + 6];
				col[4] = x[6*ii + 4 + 6];
				col[5] = x[6*ii + 5 + 6];

				//convert to eigen
				Eigen::VectorXd eigcol = Util::StdVec2Eigen(col);

				//assign
				STM.block(0,ii,6,1) = eigcol;

			} // for


			//change in STM
			Eigen::MatrixXd dSTM = jac*STM;

			//assign
			for (int ii = 0; ii < 6; ++ii) {
				Eigen::VectorXd eigcol = dSTM.block(0,ii,6,1);

				dxdt[6*ii + 0 + 6] = eigcol[0];
				dxdt[6*ii + 1 + 6] = eigcol[1];
				dxdt[6*ii + 2 + 6] = eigcol[2];
				dxdt[6*ii + 3 + 6] = eigcol[3];
				dxdt[6*ii + 4 + 6] = eigcol[4];
				dxdt[6*ii + 5 + 6] = eigcol[5];

			} // for


			
		} //fi


		//change in state
		dxdt[0] = vel[0];
		dxdt[1] = vel[1];
		dxdt[2] = vel[2];
		dxdt[3] = accel[0];
		dxdt[4] = accel[1];
		dxdt[5] = accel[2];

		// std::cout << "simple accel: \n" << accel[0] << "\n" << accel[1] << "\n" << accel[2] << "\n";

		// exit(0);

	} //operator()

	Eigen::Vector2d Propagator::GetRangeAndRate(Eigen::Vector3d pos_station_ecef, tof){

		//locals
		double JD_UTC = this->t_JD_ + this->t_/(24.0*60.0*60.0); //current UTC Julian Date in days
		double earthrot = this->earthrotationspeed_;
		Eigen::Vector3d pos_craft = this->pos_;
		Eigen::Vector3d vel_craft = this->vel_;
		std::vector<Eigen::Matrix3d> matvec(4, Eigen::Matrix3d::Identity(3,3));
		Eigen::Matrix3d R_ecef2eci = Util::ECEF2ECI(JD_UTC, this->gravmodel_->nut80ptr_, this->gravmodel_->iau1980ptr_, &matvec); //matrix rotating from ECEF2ECI


		//convert station position to ECI
		Eigen::Vector3d pos_station = R_ecef2eci*pos_station_ecef;

		//initialize measurement
		Eigen::Vector2d z;

		//station ECI velocity
		// Eigen::Vector3d vel_station;
		// vel_station[0] = -1.0*earthrot*pos_station[1];
		// vel_station[1] = earthrot*pos_station[0];
		// vel_station[2] = 0.0;

		// std::cout << "vel1: \n" << vel_station << "\n";

		// Vallado IAU-76/FK5 transformation for velocity
		//convert station position to PEF 
		Eigen::Vector3d pos_station_pef = matvec[3]*pos_station_ecef;

		//modify earth rate
		double LOD = matvec[4](0,0);
		earthrot = earthrot*(1.0 - LOD/86400.0);

		//station PEF velocity
		Eigen::Vector3d vel_station_pef;
		vel_station_pef[0] = -1.0*earthrot*pos_station_pef[1];
		vel_station_pef[1] = earthrot*pos_station_pef[0];
		vel_station_pef[2] = 0.0;

		//station ECI velocity
		Eigen::Vector3d vel_station_eci = Util::OrthoNormMat(matvec[0]*matvec[1]*matvec[2])*vel_station_pef;

		// std::cout << "vel2: \n" << R_ecef2eci*vel_station_ecef << "\n";

		// exit(0);

		//relative position and velocity
		Eigen::Vector3d rel_pos = pos_craft - pos_station + tof*vel_station_eci;
		Eigen::Vector3d rel_vel  = vel_craft - vel_station_eci;

		// std::cout << "station ECI vel: \n" << vel_station_eci << "\n";
		// exit(0);

		//range
		z[0] = rel_pos.norm();

		//range rate
		z[1] = rel_pos.dot(rel_vel)/z[0];

		return z;

	} //GetRangeAndRate

	Eigen::MatrixXd Propagator::GetRangeAndRateJac(Eigen::Vector3d pos_station_ecef){

		//locals
		double JD_UTC = this->t_JD_ + this->t_/(24.0*60.0*60.0); //current UTC Julian Date in days
		double earthrot = this->earthrotationspeed_;
		Eigen::Vector3d pos_craft = this->pos_;
		Eigen::Vector3d vel_craft = this->vel_;
		std::vector<Eigen::Matrix3d> matvec(4, Eigen::Matrix3d::Identity(3,3));
		Eigen::Matrix3d R_ecef2eci = Util::ECEF2ECI(JD_UTC, this->gravmodel_->nut80ptr_, this->gravmodel_->iau1980ptr_, &matvec); //matrix rotating from ECEF2ECI

		//convert station position to ECI
		Eigen::Vector3d pos_station = R_ecef2eci*pos_station_ecef;

		//station ECI velocity
		Eigen::Vector3d vel_station;
		vel_station[0] = -1.0*earthrot*pos_station[1];
		vel_station[1] = earthrot*pos_station[0];
		vel_station[2] = 0.0;

		//relative position and velocity
		Eigen::Vector3d rel_pos = pos_craft - pos_station;
		Eigen::Vector3d rel_vel  = vel_craft - vel_station;

		//initialize H
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2,6);

		//range partials
		H(0,0) = rel_pos[0]/rel_pos.norm();
		H(0,1) = rel_pos[1]/rel_pos.norm();
		H(0,2) = rel_pos[2]/rel_pos.norm();

		//range rate partials
		double rel32 = pow(pow(rel_pos[0],2.0) + pow(rel_pos[1],2.0) + pow(rel_pos[2],2.0), 1.5);
		H(1,0) = (rel_vel[0]*pow(rel_pos[1],2.0) - rel_pos[0]*rel_vel[1]*rel_pos[1] + rel_vel[0]*pow(rel_pos[2],2.0) - rel_pos[0]*rel_vel[2]*rel_pos[2])/rel32;
		H(1,1) = (rel_vel[1]*pow(rel_pos[0],2.0) - rel_pos[1]*rel_vel[0]*rel_pos[0] + rel_vel[1]*pow(rel_pos[2],2.0) - rel_pos[1]*rel_vel[2]*rel_pos[2])/rel32;
		H(1,2) = (rel_vel[2]*pow(rel_pos[1],2.0) - rel_pos[2]*rel_vel[0]*rel_pos[0] + rel_vel[2]*pow(rel_pos[1],2.0) - rel_pos[2]*rel_vel[1]*rel_pos[1])/rel32;
		H(1,3) = H(0,0);
		H(1,4) = H(0,1);
		H(1,5) = H(0,2);

		return H;
	} //GetRangeAndRateJac()




} //namespace VehicleState