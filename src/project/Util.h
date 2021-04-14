#ifndef OD_SRC_PROJECT_UTIL_H_
#define OD_SRC_PROJECT_UTIL_H_

#include <Eigen/Dense> //matricies and vectors
#include <vector> // basic cpp vectors
#include <string> //filenames

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

	//Create the jacobian for a gravity field with J2 and J3 and drag
	Eigen::MatrixXd GetGravJac(Eigen::Vector3d pos, Eigen::Vector3d vel, double Rearth, double wEarth, double mu, double J2, double J3, double drag_coeff);

	//convert from ECEF (ITRF) [km] to ECI (ICRF) [km] using IAU-76/FK5 and julian date (UTC)
	//	sourced from Vallado
	Eigen::Matrix3d ECEF2ECI(double JD_UTC, Eigen::MatrixXd* nut80ptr, Eigen::MatrixXd* iau1980ptr);

	//load comma seperated value data file
	Eigen::MatrixXd LoadDatFile(std::string file, int rows, int cols);

	//raise a 3x3 matrix to a given power
	Eigen::Matrix3d MatPOW(Eigen::Matrix3d mat, int n);

	//orthonomalize a 3x3 matrix
	Eigen::Matrix3d OrthoNormMat(Eigen::Matrix3d mat);

	//convert a julian date in day month hour etc into a JD [days] (vallado algo 14)
	double JulianDateNatural2JD(double year, double month, double day, double hour, double min, double sec);

	//class for EGM-96 Gravity Model
	class EGM96Grav{

	public:

		//parameters
		double mu_;
		double Rearth_; //radius of earth
		Eigen::MatrixXd C_; //un-normalized C coeffs
		Eigen::MatrixXd C_norm_; //normalized C coeffs
		Eigen::MatrixXd S_; //un-normalized S coeffs
		Eigen::MatrixXd S_norm_; //normalized S coeffs
		Eigen::MatrixXd* nut80ptr_; //nutation matrix pointer
		Eigen::MatrixXd* iau1980ptr_; //time matrix pointer

		//constructor
		EGM96Grav();

		//Load normalized coeffs from a csv file
		void LoadNormCoeffs(std::string C_file, std::string S_file);

		//convert normalized coefficients to un-normalized ones
		void NormCoeffs2Reg();

		//get acceleration
		Eigen::Vector3d GetGravAccel(Eigen::Vector3d pos, double JD_UTC);

	}; //class EGM96Grav

	//cost function for nonlinear optimization of residuals across all measurements
	double GetCost(Eigen::MatrixXd x0, Eigen::MatrixXd bias);


} //namespace Util


#endif