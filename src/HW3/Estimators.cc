#include "Estimators.h"
#include <Eigen/Dense> // vectors

namespace Estimate{

	LS::LS(){};

	void LS::CalcEstimate(Eigen::VectorXd y){

		//extract locals
		Eigen::MatrixXd H = this->H_;

		//invert R
		Eigen::MatrixXd Rinv = this->R_.inverse();

		//covariance
		this->Phat_ = (H.transpose()*Rinv*H).inverse();

		//estimate
		this->xhat_ = this->Phat_*H.transpose()*Rinv*y;


	} //Estimate


} //namespace Estimate