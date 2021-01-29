#include <iostream> //outputs
#include <boost/numeric/odeint.hpp> //integrator
#include <math.h> // trig functions
#include <Eigen/Dense> // vectors
#include "matplotlibcpp.h" // plotting

// Define the analytic solution for the harmonic oscillator
Eigen::VectorXf AnalyticOscill(float A, float ratio, Eigen::VectorXf t, float phi){

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

int main() {
	
	//initialize variables
	Eigen::VectorXf t = Eigen::VectorXf::LinSpaced(101, 0, 20);
	float A = 1.34;
	float phi = M_PI/3.0;
	float ratio = 1.0;

	//call analytic function
	Eigen::VectorXf y_analyt = AnalyticOscill(A, ratio, t, phi);

	std::cout << y_analyt << "\n";

}