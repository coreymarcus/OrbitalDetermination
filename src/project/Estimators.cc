#include "Estimators.h"
#include <Eigen/Dense> // vectors
#include <iostream> //outputs

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

	EKF::EKF(){};

	void EKF::CalcEstimate(Eigen::VectorXd y){

		//extract locals
		Eigen::MatrixXd H = this->H_;
		Eigen::MatrixXd R = this->R_;
		Eigen::MatrixXd Pbar = this->Phat_;

		// find kalman gain
		Eigen::MatrixXd Pxy = Pbar*H.transpose();
		Eigen::MatrixXd Pyy = H*Pbar*H.transpose() + R;
		Eigen::MatrixXd K = Pxy*Pyy.inverse();

		//covariance
		Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(this->n_,this->n_) - K*H;
		this->Phat_ = Z*Pbar*Z.transpose() + K*R*K.transpose();

		//estimate
		this->xhat_ = K*y;


	} //Estimate

	UKF::UKF(){};

	void UKF::GetSigmaPoints(){

		//Extract Locals
		Eigen::MatrixXd Phat = this->Phat_;
		Eigen::VectorXd xhat = this->xhat_;
		int n = this->n_;
		double npk = (double) n + this->k_;

		//initialize Xi and weights
		Eigen::MatrixXd Xi = Eigen::MatrixXd::Zero(n, 2*n + 1);
		Eigen::VectorXd w = Eigen::VectorXd::Zero(2*n + 1);

		

		//get the matrix square root
		Eigen::MatrixXd S = Phat.llt().matrixL();

		//get non-central points
		for (int ii = 0; ii < n; ++ii) {
			Xi.block(0,ii,n,1) = xhat - sqrt(npk)*S.block(0,ii,n,1);
			Xi.block(0,ii + n,n,1) = xhat + sqrt(npk)*S.block(0,ii,n,1);
			w(ii) = 0.5/npk;
			w(ii+n) = 0.5/npk;
			
		}

		//central point
		Xi.block(0,2*n,n,1) = xhat;
		w(2*n) = this->k_ / npk;

		//assign
		this->Xi_ = Xi;
		this->w_ = w;


	} //PropEstimate

	void UKF::SigmaPts2Estimate(){

		//extract locals
		Eigen::MatrixXd Xi = this->Xi_;
		Eigen::VectorXd w = this->w_;
		int n = this->n_;

		//initialize variables
		Eigen::VectorXd xhat = Eigen::VectorXd::Zero(n);
		Eigen::MatrixXd Phat = Eigen::MatrixXd::Zero(n,n);

		//get estimate
		for (int i = 0; i < (2*n + 1); ++i) {
			xhat += w(i)*Xi.block(0,i,n,1);
		}

		//get covariance
		for (int i = 0; i < (2*n + 1); ++i) {
			Eigen::MatrixXd d = Xi.block(0,i,n,1) - xhat;
			Phat += w(i)*d*d.transpose();
		}

		//assign
		this->Phat_ = Phat;
		this->xhat_ = xhat;


	} //SigmaPts2Estimate()

	void UKF::CalcEstimate(Eigen::VectorXd y){

	} //CalcEstimate


} //namespace Estimate