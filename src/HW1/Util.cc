#include "Util.h"
#include <iostream> //outputs
#include <Eigen/Dense> // vectors
#include <math.h> // trig functions
#include <vector> // basic cpp vectors

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

} //namespace Util