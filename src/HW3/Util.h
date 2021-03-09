#ifndef OD_SRC_HW1_UTIL_H_
#define OD_SRC_HW1_UTIL_H_

#include <Eigen/Dense> //matricies and vectors
#include <vector> // basic cpp vectors

namespace Util {

	/* Create a 3x3 rotation matrix specified by a series of angles and an order.
		Rotation is a frame rotation
		angles specified like "123" or "313"
	*/
	Eigen::Matrix3d Angle2RotM(Eigen::Vector3d angles, std::vector<int> order);

	//quickly convert from eigen vector to std::vector
	std::vector<double> EigenVec2Std(Eigen::VectorXd vec1);

	//quickly convert from std::vector to eigen vector
	Eigen::VectorXd StdVec2Eigen(std::vector<double> vec1);

	//write eigen matrix or vector to csv file
	void Eigen2csv(std::string file, Eigen::MatrixXd mat);

	//Create the jacobian for a simple newtonian gravity field
	Eigen::MatrixXd GetGravJac(Eigen::Vector3d pos, double mu);

} //namespace Util


#endif