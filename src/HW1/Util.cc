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
			vec2[ii] = vec1(ii);
		}

		//output
		return vec2;
	} //EigenVec2Std

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

} //namespace Util