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

	KF::KF(){};

	void KF::CalcEstimate(Eigen::VectorXd y){

		//extract locals
		Eigen::MatrixXd H = this->H_;
		Eigen::MatrixXd R = this->R_;
		Eigen::MatrixXd Pbar = this->Pbar_;
		int N_x = this->mu_.size();

		// find kalman gain
		Eigen::MatrixXd Pxy = Pbar*H.transpose();
		Eigen::MatrixXd Pyy = H*Pbar*H.transpose() + R;
		Eigen::MatrixXd K = Pxy*Pyy.inverse();

		//covariance
		Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(N_x,N_x) - K*H;
		this->Phat_ = Z*Pbar*Z.transpose() + K*R*K.transpose();

		//estimate
		this->xhat_ = this->mu_ + K*(y - H*this->mu_);


	} //Estimate


} //namespace Estimate