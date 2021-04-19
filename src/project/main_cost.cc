#include <Eigen/Dense> // vectors
#include "VehicleState.h" // headerfile containing propagator class
#include "Util.h" // utility functions



int main() {

	//load data
	Eigen::MatrixXd xeval = Util::LoadDatFile("../data/xeval_project.csv", 13, 1);

	double test = Util::GetCost(xeval.block(0,0,6,1), xeval.block(6,0,6,1), xeval(12,0));
	Eigen::MatrixXd testmat = Eigen::MatrixXd::Zero(1,1);
	testmat(0,0) = test;

	Util::Eigen2csv("../data/evalcost_project.csv", testmat);

	return 0;
}