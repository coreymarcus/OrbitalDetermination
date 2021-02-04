#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "matplotlibcpp.h" // plotting

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

//namespaces, don't ever use name space, simply call them aliases
namespace plt = matplotlibcpp;
namespace ode = boost::numeric::odeint;

// Define the analytic solution for the harmonic oscillator
Eigen::VectorXf AnalyticOscill(double A, double ratio, Eigen::VectorXf t, double phi){

	//length of t
	int n = t.size();

	//initialize output
	Eigen::VectorXf y(n);

	//loop
	for (int ii = 0; ii < n; ++ii) {
		y(ii) = A*cos( sqrt(ratio)*t(ii) + phi );
	}

	//output
	return  y;
}

//quickly convert from eigen vector to std:vector
std::vector<double> EigenVec2Std(Eigen::VectorXf vec1){
	
	//locals
	int n = vec1.size();
	std::vector<double> vec2(n,0);

	//loop
	for (int ii = 0; ii < n; ++ii){
		vec2[ii] = vec1(ii);
	}

	//output
	return vec2;
}

// harmonic oscillator class
class HarmOscillator {

public:

	//Data Members
	double kmratio;

	//Member Functions
	void operator() ( const state_type &x , state_type &dxdt , const double /* t */ ) {

		dxdt[0] = x[1];
		dxdt[1] = -kmratio*x[0];
	}

};

int main() {
	
	//initialize variables
	int n = 101; //number of points
	Eigen::VectorXf t = Eigen::VectorXf::LinSpaced(n, 0, 20);
	double A = 1.34;
	double phi = M_PI/3.0;
	double kmratio = 1.0;

	//initialize ode stuff
	// this is the ode solver we would like to use
	typedef ode::runge_kutta_dopri5< state_type > error_stepper_type;
	//this is used to control error parameters
	typedef ode::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
	HarmOscillator fun;
	fun.kmratio = kmratio; // initialize our dynamics
	const double dt_var = 0.1; /* this is the initial guess for a variable step
		size for the ode solver */
	const double reltol = 1*pow(10,-12); //relative tolerance
	const double abstol = 1*pow(10,-20); //absolute tolerance

	//initial conditions
	state_type x0(2);
	x0[0] = A*cos(phi);
	x0[1] = -A*sqrt(kmratio)*sin(phi);

	//call analytic function
	Eigen::VectorXf y_analyt = AnalyticOscill(A, kmratio, t, phi);

	//integrate numerically
	state_type xint = x0; //state for integration
	std::vector<double> y_num(n,0); //storage for numerically integrated state
	y_num[0] = x0[0];
	for (int ii = 1; ii < n; ++ii)	{

		//get times
		double t1 = t[ii-1];
		double t2 = t[ii];

		//integrate
		ode::integrate_adaptive( ode::make_controlled( abstol , reltol ,
			error_stepper_type() ) , fun , xint , t1 , t2 , dt_var );

		//extract result
		y_num[ii] = xint[0];
		
	}

	//find error
	std::vector<double> interr(n,0);
	for (int ii = 0; ii < n; ++ii){
		interr[ii] = y_num[ii] - y_analyt[ii];
	}

	//plot results
	std::vector<double> t_std = EigenVec2Std(t);
	std::vector<double> y_analyt_std = EigenVec2Std(y_analyt);

	plt::figure(1);
	plt::named_plot("Analytic",t_std, y_analyt_std);
	plt::title("Homework 1 - Plot 1");
	plt::xlabel("Time (s)");
	plt::ylabel("Analytic System Response");
	// plt::legend();
    

    plt::figure(2);
    plt::named_plot("Error",t_std, interr);
	plt::title("Homework 1 - Plot 2");
	plt::xlabel("Time (s)");
	plt::ylabel("Integration Error");

	plt::show();

}