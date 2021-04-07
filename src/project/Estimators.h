#ifndef OD_SRC_PROJECT_ESTIMATORS_H_
#define OD_SRC_PROJECT_ESTIMATORS_H_

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

	//Extended Kalman Filter for x(k+1) = f( x(k) ) + v; y(k) = h( x(k) ) + w
	//estimates deviations from nominal state!!!
	class EKF{

	public:

		//parameterss
		Eigen::MatrixXd F_; //jacobian of f()
		Eigen::MatrixXd H_; //jacobian of h()
		Eigen::MatrixXd Q_; //process noise covariance
		Eigen::MatrixXd R_; //measurement noise covariance
		int n_;

		//Estimates
		Eigen::MatrixXd Phat_; //estimate covariance
		Eigen::VectorXd xhat_;

		// Constructor
		EKF();

		// performs an estimate given measurement y of deviation (y - h(xhat))
		void CalcEstimate(Eigen::VectorXd y);

	}; //class KF

	// Unscented Kalman Filter
	class UKF{

	public:

		//parameters
		double k_; //unnormalized weight of the central point for '2n+1' sigma points

		Eigen::MatrixXd Q_; //process noise covariance
		Eigen::MatrixXd R_; //measurement noise covariance
		int n_; //size of the state
		Eigen::MatrixXd Xi_; //storage for sigma points
		Eigen::VectorXd w_; //storage for weights

		//estimates
		Eigen::MatrixXd Phat_;
		Eigen::VectorXd xhat_;

		// Constructor
		UKF();

		// populate sigma points for the current estimate
		void GetSigmaPoints();

		// use the sigma points to calculate the new estimate
		void SigmaPts2Estimate();

		// Performs an update given measurement  y
		void CalcEstimate(Eigen::VectorXd y);

	};

} //namespace Estimate

#endif