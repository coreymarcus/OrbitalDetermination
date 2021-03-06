#include "Util.h"
#include <iostream> //outputs
#include <Eigen/Dense> // vectors
#include <math.h> // trig functions
#include <vector> // basic cpp vectors
#include <fstream> // writing to files
#include <iomanip>  // std::setprecision()

namespace Util {

	Eigen::Matrix3d Angle2RotM(Eigen::Vector3d angles, std::vector<int> order){

		//initialize matrix
		Eigen::Matrix3d R = Eigen::Matrix3d::Identity(3,3);

		for (int ii = 0; ii < 3; ++ii){
			
			//iteration matrix
			Eigen::Matrix3d R_iter = Eigen::Matrix3d::Identity(3,3);

			//switch case
			switch(order[ii]) {
				case 1 : {
					R_iter(1,1) = cos(angles[ii]);
					R_iter(1,2) = sin(angles[ii]);
					R_iter(2,1) = -sin(angles[ii]);
					R_iter(2,2) = cos(angles[ii]);
					break;
				} case 2 : {
					R_iter(0,0) = cos(angles[ii]);
					R_iter(0,2) = -sin(angles[ii]);
					R_iter(2,0) = sin(angles[ii]);
					R_iter(2,2) = cos(angles[ii]);
					break;
				} case 3 : {
					R_iter(0,0) = cos(angles[ii]);
					R_iter(0,1) = sin(angles[ii]);
					R_iter(1,0) = -sin(angles[ii]);
					R_iter(1,1) = cos(angles[ii]);
					break;
				} default : {
					std::cout << "Error: Bad Angle Selection!" << std::endl;
				}
			} // switch

			//multiply
			R = R_iter*R;

		} // for

		//output
		return R.transpose(); // this is a transpose to make it frame rot instead of vector rot

	} // Angle2RotM

	std::vector<double> EigenVec2Std(Eigen::VectorXd vec1){
		
		//locals
		int n = vec1.size();
		std::vector<double> vec2(n,0);

		//loop
		for (int ii = 0; ii < n; ++ii){
			vec2[ii] = vec1[ii];
		}

		//output
		return vec2;
	} //EigenVec2Std

	Eigen::VectorXd StdVec2Eigen(std::vector<double> vec1){

		//locals
		int n = vec1.size();
		Eigen::VectorXd vec2(n);

		//loop
		for (int ii = 0; ii < n; ++ii){
			vec2[ii] = vec1[ii];
		}

		//output
		return vec2;

	} //StdVec2Eigen

	void Eigen2csv(std::string file, Eigen::MatrixXd mat){

		//locals
		int height = mat.rows();
		int width = mat.cols();

		//open file
		std::ofstream f;
		f.open(file);
		if(!f.is_open()){
			std::cout << "Error: Unable to open file!" << std::endl;
		}

		//set precision
		f << std::fixed << std::setprecision(17);

		//write
		for (int ii = 0; ii < height; ++ii) {
			for (int jj = 0; jj < width; ++jj) {
			
				// write element
				f << mat(ii,jj);

				//write comma
				if(jj != (width-1)){
					f << ", ";
				}

			}

			//newline
			f << std::endl;
		}

		f.close();


	} //Eigen2csv

	Eigen::MatrixXd GetGravJac(Eigen::Vector3d pos, double mu){

		//extract locals
		double x = pos[0];
		double y = pos[1];
		double z = pos[2];

		//initialize matrix
		Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(6,6);

		//velocity component
		jac.block(0,3,3,3) = Eigen::Matrix3d::Identity(3,3);

		//helper functions
		double fun4 = pow( x*x + y*y + z*z , 2.5);
		double fun1 = 3*mu*y*z/fun4;
		double fun2 = 3*mu*x*z/fun4;
		double fun3 = 3*mu*x*y/fun4;

		//assign
		jac(3,0) = -mu*(-2*x*x + y*y + z*z)/fun4;
		jac(3,1) = fun3;
		jac(3,2) = fun2;

		jac(4,0) = fun3;
		jac(4,1) = -mu*(x*x - 2*y*y + z*z)/fun4;
		jac(4,2) = fun1;

		jac(5,0) = fun2;
		jac(5,1) = fun1;
		jac(5,3) = -mu*(x*x + y*y - 2*z*z)/fun4;

		//return
		return jac;

	} //GetGravJac

} //namespace Util