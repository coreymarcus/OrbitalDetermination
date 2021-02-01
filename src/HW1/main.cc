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

//quickly convert from eigen vector to std:vector
std::vector<float> EigenVec2Std(Eigen::VectorXf vec1){
	
	//locals
	int n = vec1.size();
	std::vector<float> vec2(n,0);

	//loop
	for (int ii = 0; ii < n; ++ii){
		vec2[ii] = vec1(ii);
	}

	//output
	return vec2;
}

//namespaces, don't ever use name space, simply call them aliases
namespace plt = matplotlibcpp;

int main() {
	
	//initialize variables
	Eigen::VectorXf t = Eigen::VectorXf::LinSpaced(101, 0, 20);
	float A = 1.34;
	float phi = M_PI/3.0;
	float ratio = 1.0;

	//call analytic function
	Eigen::VectorXf y_analyt = AnalyticOscill(A, ratio, t, phi);


	//plot results
	std::vector<float> t_std = EigenVec2Std(t);
	std::vector<float> y_analyt_std = EigenVec2Std(y_analyt);
	std::vector<float> y_analyt_std2 = EigenVec2Std(2*y_analyt);

	plt::named_plot("Analytic",t_std, y_analyt_std);
	plt::named_plot("test2",t_std, y_analyt_std2);
	plt::title("Homework 1 - Plot 1");
	plt::xlabel("Time (s)");
	plt::ylabel("System Response");
	plt::legend();
    plt::show();

}