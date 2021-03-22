#ifndef OD_SRC_HW4_ESTIMATORS_H_
#define OD_SRC_HW4_ESTIMATORS_H_

#include <Eigen/Dense>

namespace Estimate{

	//linear least squares estimator for y = Hx + w
	class LS{

	public:

		//parameters
		Eigen::MatrixXd H_; //measurement mapping
		Eigen::MatrixXd R_; //measurement noise covariance
		Eigen::MatrixXd Pbar_; //prior covariance of x
		Eigen::VectorXd mu_; //prior mean of x

		//Estimates
		Eigen::VectorXd xhat_; //estimate of x
		Eigen::MatrixXd Phat_; //estimate covariance

		// Constructor
		LS();

		// performs an estimate given measurement y
		void CalcEstimate(Eigen::VectorXd y);

	}; //class LS

	//Kalman Filter for y = Hx + w (static system!)
	class KF{

	public:

		//parameters
		Eigen::MatrixXd H_; //measurement mapping
		Eigen::MatrixXd R_; //measurement noise covariance
		Eigen::MatrixXd Pbar_; //prior covariance of x
		Eigen::VectorXd mu_; //prior mean of x

		//Estimates
		Eigen::VectorXd xhat_; //estimate of x
		Eigen::MatrixXd Phat_; //estimate covariance

		// Constructor
		KF();

		// performs an estimate given measurement y
		void CalcEstimate(Eigen::VectorXd y);

	}; //class KF

} //namespace Estimate

#endif