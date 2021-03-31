#include "Util.h"
#include <iostream> //outputs
#include <Eigen/Dense> // vectors
#include <math.h> // trig functions
#include <vector> // basic cpp vectors
#include <fstream> // writing to files
#include <iomanip>  // std::setprecision()
#include <string> //filenames
#include <boost/math/special_functions/legendre.hpp> //legendre polynomials for gravity

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
		return OrthoNormMat(R); // not super sure if this should actually be R.transpose()

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
		jac(5,2) = -mu*(x*x + y*y - 2*z*z)/fun4;

		//return
		return jac;

	} //GetGravJac

	Eigen::Matrix3d ECEF2ECI(double JD_UTC, Eigen::MatrixXd* nut80ptr, Eigen::MatrixXd* iau1980ptr){

		//********* deal with time ******************

		//J2000 time
		double J2000 = 2451545.0;

		//offset to get mean julian dates
		double MJD_offset = 2400000.5;

		//mean julian date, J2000
		double MJD_J2000 = J2000 - MJD_offset;

		//current mean JD
		double MJD_UTC = JD_UTC - MJD_offset;

		double MJD_TAI = MJD_UTC + 37.0/(3600.0*24.0);
		double MJD_TT = MJD_TAI + 32.184/(3600.0*24.0);

		// find the change in time since J2000.0 (T_TT)
		double tsinceJ2000 = (MJD_TT - MJD_J2000)/36525.0;

		//find the remainder in the current JD
		double JD_UTC_rem = MJD_UTC - floor(MJD_UTC);

		// std::cout << "MJD_J2000: " << MJD_J2000 << "\n";
		// std::cout << "MJD_UTC: " << MJD_UTC << "\n";
		// std::cout << "MJD_TAI: " << MJD_TAI << "\n";
		// std::cout << "MJD_TT: " << MJD_TT << "\n";
		// std::cout << "JD_UTC_rem: " << JD_UTC_rem << "\n";

		//we will manually find values from iau1980.txt because lazy
		// std::cout << "You need to manually find this in iau1980.txt \n";
		// std::cout << "floor(MJD_UTC): " << floor(MJD_UTC) << "\n";

		//manual search found this:
		// 1712 1 58088.00 I  0.124136 0.000022  0.236728 0.000035  I 0.2485001 0.0000084  1.5396 0.0061  I  -104.226    0.322    -8.561    0.160  0.124126  0.236687  0.2485227  -104.524    -8.685  
		// 1712 2 58089.00 I  0.121708 0.000026  0.236266 0.000028  I 0.2469699 0.0000084  1.4837 0.0058  I  -104.143    0.322    -8.430    0.160  0.121669  0.236396  0.2469833  -105.014    -8.516
		// double PMx1 = 0.124126;
		// double PMx2 = 0.121669;
		// double PMy1 = 0.236687;
		// double PMy2 = 0.236396;
		// double UT1_UTC1 = 0.24852270;
		// double UT1_UTC2 = 0.24698330;

		//get initial index of desired time
		int tidx = round(floor(MJD_UTC) - iau1980ptr->coeff(0,0));

		// std::cout << JD_UTC << "\n";

		//extract values
		double PMx1 = iau1980ptr->coeff(tidx,1);
		double PMx2 = iau1980ptr->coeff(tidx+1,1);
		double PMy1 = iau1980ptr->coeff(tidx,2);
		double PMy2 = iau1980ptr->coeff(tidx+1,2);
		double UT1_UTC1 = iau1980ptr->coeff(tidx,3);
		double UT1_UTC2 = iau1980ptr->coeff(tidx+1,3);

		//******** find rotation from precession ******************

		//now, find the three angles in degrees
		double zeta_deg = 2306.2181*tsinceJ2000 + 0.30188*pow(tsinceJ2000,2) + 0.017998*pow(tsinceJ2000,3); 
		double theta_deg = 2004.3109*tsinceJ2000 - 0.42665*pow(tsinceJ2000,2) - 0.041833*pow(tsinceJ2000,3); 
		double z_deg = 2306.2181*tsinceJ2000 + 1.09468*pow(tsinceJ2000,2) + 0.018203*pow(tsinceJ2000,3);

		//convert angles to radians
		double zeta = zeta_deg*M_PI/(3600.0*180.0);
		double theta = theta_deg*M_PI/(3600.0*180.0);
		double z = z_deg*M_PI/(3600.0*180.0);

		// std::cout << "tsinceJ2000: " << tsinceJ2000 << "\n";
		// std::cout << "zeta: " << zeta << "\n";
		// std::cout << "theta: " << theta << "\n";
		// std::cout << "z: " << z << "\n";

		//create matrix
		Eigen::Vector3d angles(z,-theta,zeta);
		std::vector<int> order = {3,2,3};
		Eigen::Matrix3d P = Angle2RotM(angles, order);

		//******** find rotation from nutation ******************

		// load the parameters
		// Eigen::MatrixXd nut80 = LoadDatFile("../data/nut80.csv", 106, 10);

		//find epsbar1980
		double epsbar1980 = (84381.448 - 46.8150*tsinceJ2000 - 0.00059*pow(tsinceJ2000,2) + 0.001813*pow(tsinceJ2000,3))*M_PI/(3600.0*180.0);

		//find values from 3-82 (corrected according to errata)
		double r = 360.0;
		double M_moon = 134.96298139 + (1325*r + 198.8673981)*tsinceJ2000 + 0.0086972*pow(tsinceJ2000,2) + 1.78*pow(10.0,-5)*pow(tsinceJ2000,3);
		double M_cirle = 357.52772333 + (99*r + 359.0503400)*tsinceJ2000 - 0.0001603*pow(tsinceJ2000,2) - 3.3*pow(10.0,-6)*pow(tsinceJ2000,3);
		double uM_moon =  93.27191028 + (1342*r + 82.0175381)*tsinceJ2000 - 0.0036825*pow(tsinceJ2000,2) + 3.1*pow(10.0,-6)*pow(tsinceJ2000,3);
		double D_circle = 297.85036306 + (1236*r + 307.1114800)*tsinceJ2000 - 0.0019142*pow(tsinceJ2000,2) + 5.3*pow(10.0,-6)*pow(tsinceJ2000,3);
		double Ohm_moon = 125.04452222 - (5*r + 134.1362608)*tsinceJ2000 +  0.0020708*pow(tsinceJ2000,2) + 2.2*pow(10.0,-6)*pow(tsinceJ2000,3);

		//convert to radians
		M_moon = M_moon*M_PI/180.0;
		M_cirle = M_cirle*M_PI/180.0;
		uM_moon = uM_moon*M_PI/180.0;
		D_circle = D_circle*M_PI/180.0;
		Ohm_moon = Ohm_moon*M_PI/180.0;

		// std::cout << "M_moon: " << M_moon << "\n";
		// std::cout << "M_cirle: " << M_cirle << "\n";
		// std::cout << "uM_moon: " << uM_moon << "\n";
		// std::cout << "D_circle: " << D_circle << "\n";
		// std::cout << "Ohm_moon: " << Ohm_moon << "\n";

		//sum
		double deltaPsi_1980 = 0.0;
		double deltaEsp_1980 = 0.0;
		for (int ii = 0; ii < 106; ++ii) {

			//extract locals
			double an1 = nut80ptr->coeff(ii,0);
			double an2 = nut80ptr->coeff(ii,1);
			double an3 = nut80ptr->coeff(ii,2);
			double an4 = nut80ptr->coeff(ii,3);
			double an5 = nut80ptr->coeff(ii,4);
			double A = nut80ptr->coeff(ii,5);
			double B = nut80ptr->coeff(ii,6);
			double C = nut80ptr->coeff(ii,7);
			double D = nut80ptr->coeff(ii,8);
			double idx = nut80ptr->coeff(ii,9);

			//aPi
			double aPi = an1*M_moon + an2*M_cirle + an3*uM_moon + an4*D_circle + an5*Ohm_moon;

			//sums
			deltaPsi_1980 = deltaPsi_1980 + (A + B*tsinceJ2000)*sin(aPi);
			deltaEsp_1980 = deltaEsp_1980 + (C + D*tsinceJ2000)*cos(aPi);
		}

		//convert sums to radians
		deltaPsi_1980 = deltaPsi_1980*pow(10.0,-4)*M_PI/(3600.0*180.0);
		deltaEsp_1980 = deltaEsp_1980*pow(10.0,-4)*M_PI/(3600.0*180.0);

		// std::cout << "deltaPsi_1980: " << deltaPsi_1980 << "\n";
		// std::cout << "deltaEsp_1980: " << deltaEsp_1980 << "\n";

		//rotation
		order = {1,3,1};
		angles[0] = epsbar1980 + deltaEsp_1980;
		angles[1] = deltaPsi_1980;
		angles[2] = -epsbar1980;
		Eigen::Matrix3d N = Angle2RotM(angles,order);

		// std::cout << "angles: " << angles << "\n";

		//******** find rotation from Polar Motion ******************

		//interpolate these two values based on iau1980.txt
		double xp = PMx1 + JD_UTC_rem*(PMx2 - PMx1);
		double yp = PMy1 + JD_UTC_rem*(PMy2 - PMy1);

		//convert to radians
		xp = xp*M_PI/(3600.0*180.0);
		yp = yp*M_PI/(3600.0*180.0);

		//rotation
		order = {2,1,3};
		angles[0] = xp;
		angles[1] = yp;
		angles[2] = 0.0;
		Eigen::Matrix3d W = Angle2RotM(angles,order);

		//******** find rotation from Earth Rotation ******************

		//interpolate change in UTC
		double dUTC = UT1_UTC1 + JD_UTC_rem*(UT1_UTC2 - UT1_UTC1);

		//correct modified julian date
		double MJD_UT1 = MJD_UTC + dUTC/(24.0*3600.0);

		//change in time
		double dT = MJD_UT1 - MJD_J2000;

		//Greenwich mean sidereal time (from tapley shuctz born)
		double theta_GMST = 4.894961212823058751375704430 + dT*(6.300388098984893552276513720 + dT*(5.075209994113591478053805523*pow(10.0,-15) - 9.253097568194335640067190688*pow(10.0,-24)*dT));

		//now we get some other crazy value
		double Eqequinox1982 = deltaPsi_1980*cos(epsbar1980) + (0.00264*sin(Ohm_moon) + 0.000063*sin(2*Ohm_moon))*M_PI/(3600.0*180.0);

		//find theta
		double theta_GAST = theta_GMST + Eqequinox1982;

		//rotation
		order = {3,1,3};
		angles[0] = -theta_GAST;
		angles[1] = 0.0;
		angles[2] = 0.0;
		Eigen::Matrix3d S = Angle2RotM(angles,order);

		// std::cout << "JD_UTC_rem: " << JD_UTC_rem << "\n";
		// std::cout << "dUTC: " << dUTC << "\n";
		// std::cout << "UT1_UTC1: " << UT1_UTC1 << "\n";
		// std::cout << "UT1_UTC2: " << UT1_UTC2 << "\n";

		//enforce orthonormalization
		P = OrthoNormMat(P);
		N = OrthoNormMat(N);
		S = OrthoNormMat(S);
		W = OrthoNormMat(W);

		//write out all the matricies
		// std::cout << "P:\n" << P << "\n";
		// std::cout << "N:\n" << N << "\n";
		// std::cout << "S:\n" << S << "\n";
		// std::cout << "W:\n" << W << "\n";

		//Matrix
		Eigen::Matrix3d Rotm = P*N*S*W;
		Rotm = OrthoNormMat(Rotm);

		//output
		return Rotm;

	} // ECEF2ECI()

	Eigen::MatrixXd LoadDatFile(std::string file, int rows, int cols){

		//initialize the matrix
		Eigen::MatrixXd dat = Eigen::MatrixXd::Zero(rows, cols);

		//create input filestream
		std::ifstream myFile(file);

		// Make sure the file is open
		if(!myFile.is_open()) throw std::runtime_error("Could not open file");

		//line and word string
		std::string line;
		std::string word;
		std::string::size_type sz;     // alias of size_t

		// cycle
	    for (int rowidx = 0; rowidx < rows; ++rowidx) {

	    	//get line and create string stream
	    	std::getline(myFile, line);
	    	std::stringstream ss(line);
	    	
	    	for (int colidx = 0; colidx < cols; ++colidx) {

	    		//extract
	    		std::getline(ss, word, ',');

	    		dat(rowidx,colidx) = std::stod(word, &sz);
	    		
	    	}
	    }

	    //output
	    return dat;


	} //LoadDatFile

	Eigen::Matrix3d matPOW(Eigen::Matrix3d mat, int n){

		//initialize
		Eigen::Matrix3d out = mat;

		//cycle
		for (int i = 1; i < n; ++i)	{
			out = out*mat;
		}

		return out;

	} //matPOW

	Eigen::Matrix3d OrthoNormMat(Eigen::Matrix3d mat){

		//extract columns
		Eigen::Vector3d x = mat.block(0,0,3,1);
		Eigen::Vector3d y = mat.block(0,1,3,1);
		Eigen::Vector3d z = mat.block(0,2,3,1);

		//error
		double err = x.dot(y);

		//orthogonal vectors
		Eigen::Vector3d x_ort = x - 0.5*err*y;
		Eigen::Vector3d y_ort = y - 0.5*err*z;
		Eigen::Vector3d z_ort = x_ort.cross(y_ort);

		//normalize
		x_ort = x_ort/x_ort.norm();
		y_ort = y_ort/y_ort.norm();
		z_ort = z_ort/z_ort.norm();

		//output
		Eigen::Matrix3d out;
		out.block(0,0,3,1) = x_ort;
		out.block(0,1,3,1) = y_ort;
		out.block(0,2,3,1) = z_ort;

		return out;


	}

	double JulianDateNatural2JD(double year, double month, double day, double hour, double min, double sec){

		double term1 = 367.0*year;
		double term2 = floor((month + 9.0)/12.0);
		double term3 = floor(275.0*month/9.0);
		double term4 = ((sec/60.0 + min)/60.0 + hour)/24.0;

		double JD = term1 - floor(7.0*(year + term2)/4.0) + term3 + day + 1721013.5 + term4;

		return JD;
	}

	EGM96Grav::EGM96Grav(){};

	void EGM96Grav::LoadNormCoeffs(std::string C_file, std::string S_file){

		this->C_norm_ = LoadDatFile(C_file, 361, 361);
		this->S_norm_ = LoadDatFile(S_file, 361, 361);

	} //LoadNormCoeffs

	void EGM96Grav::NormCoeffs2Reg(){

		//initialize matricies
		this->C_ = Eigen::MatrixXd::Zero(361,361);
		this->S_ = Eigen::MatrixXd::Zero(361,361);

		for (double ll = 2; ll < 362; ++ll) {

			for (double mm = 0; mm < ll+1; ++mm) {

				
				double k = 2.0;
				if( (int) mm == 0) {
					k = 1.0;
				}

				//intermediate values
				//Note: tgamma(i+1) = factorial(i)
				double a = tgamma(ll + mm + 1.0);
				double b = tgamma(ll - mm + 1.0);
				double c = (2.0*ll + 1.0);

				//find factor
				double Pi_lm = sqrt(a/(b*k*c));

				//update
				this->C_( (int) ll, (int) mm) = this->C_norm_( (int) ll, (int) mm)/Pi_lm;
				this->S_( (int) ll, (int) mm) = this->S_norm_( (int) ll, (int) mm)/Pi_lm;
			}
		}

	} //NormCoeffs2Reg

	Eigen::Vector3d EGM96Grav::GetGravAccel(Eigen::Vector3d pos, double JD_UTC){

		//make short namespace
		namespace bm = boost::math;

		//initialize output
		Eigen::Vector3d accel;

		//locals
		double r = pos.norm();
		double r2 = pow(r,2.0);
		double r3 = pow(r,3.0);
		double R = this->Rearth_;
		double mu = this->mu_;
		double pUpr = 0.0; //partial potential wrt radius
		double pUpphi = 0.0; //partial potential wrt satellite latitude
		double pUplamb = 0.0; //partial potential wrt satellite latitude
		double nr_ij = sqrt(pow(pos[0],2.0) + pow(pos[1],2.0));

		//matrix rotating from ECEF2ECI
		Eigen::Matrix3d R_ecef2eci = ECEF2ECI(JD_UTC, this->nut80ptr_, this->iau1980ptr_);

		//ecef position
		Eigen::Vector3d pos_ecef = R_ecef2eci.transpose()*pos;

		// std::cout << "pos_ecef: \n" << pos_ecef << "\n";

		//lat and long
		double phi = asin(pos_ecef[2]/r);
		double lamb = atan2(pos_ecef[1],pos_ecef[0]);

		// std::cout << "r: " << r << "\n phi: " << phi << "\n lamb: " << lamb << "\n";

		//fixed order gravity model but have capability to do more
		for (double ll = 2; ll < 21; ++ll) {

			for (double mm = 0; mm < ll+1; ++mm) {

		// for (double ll = 2; ll < 3; ++ll) {

		// 	for (double mm = 0; mm < 1; ++mm) {

				//locals
				double powR = pow(R/r,ll);
				double C = this->C_( (int) ll, (int) mm);
				double S = this->S_( (int) ll, (int) mm);
				double Plm = pow(-1.0,mm)*bm::legendre_p( (int) ll, (int) mm, sin(phi)); //pow(-1,m) to account for cordon shortley
				// double Plm_plus = pow(-1.0,mm+1.0)*bm::legendre_p( (int) ll, ((int) mm) + 1, sin(phi));
				double Plm_plus = pow(-1.0,mm)*bm::legendre_p( (int) ll - 1, (int) mm, sin(phi));

				// std::cout << "l: " << ll << " m: " << mm << "\n";

				pUpr += powR*(ll + 1.0)*Plm*(C*cos(mm*lamb) + S*sin(mm*lamb));
				// pUpphi += powR*(Plm_plus - mm*tan(phi)*Plm) * (C*cos(mm*lamb) + S*sin(mm*lamb));
				pUpphi += powR*(-ll*sin(phi)*(1/cos(phi))*Plm + (ll+mm)*(1/cos(phi))*Plm_plus) * (C*cos(mm*lamb) + S*sin(mm*lamb));
				pUplamb += powR*mm*Plm*(S*cos(mm*lamb) - C*sin(mm*lamb));

			}

		}

		//additional factors outside of sum for partials
		// pUpr = pUpr*(-mu/pow(r,2.0));
		// pUpphi = pUpphi*mu/r;
		// pUplamb = pUplamb*mu/r;
		pUpr = pUpr*(-mu/pow(r,2.0));
		pUpphi = pUpphi*mu/pow(r,2.0);
		pUplamb = pUplamb*mu/pow(r,2.0)/cos(phi);

		// pUpphi = -7.0036240305764034*pow(10.0,-8.0);

		// std::cout << "pUpr: " << pUpr << "\n" << "pUpphi: " << pUpphi << "\n" << "pUplamb: " << pUplamb << "\n";

		//calculate components of acceleration
		// accel[0] = (pUpr/r - pos[2]*pUpphi/(r2*nr_ij))*pos[0] - pUplamb*pos[1]/(pow(pos[0],2.0) + pow(pos[1],2.0)) - mu*pos[0]/r3;
		// accel[1] = (pUpr/r - pos[2]*pUpphi/(r2*nr_ij))*pos[1] + pUplamb*pos[0]/(pow(pos[0],2.0) + pow(pos[1],2.0)) - mu*pos[1]/r3;
		// accel[2] = pUpr*pos[2]/r + nr_ij*pUpphi/r2 - mu*pos[2]/r3;
		accel[0] = pUpr;
		accel[1] = pUplamb;
		accel[2] = pUpphi;

		//rotation from spherical to cartesian
		Eigen::Vector3d angles = {phi, -lamb, 0.0};
		std::vector<int> order = {2,3,1};
		Eigen::Matrix3d Rsp2cart = Angle2RotM(angles, order);

		//we also need the central acceleration in ECI
		Eigen::Vector3d accel_central = -(mu/r3)*pos;

		//rotate accel into ECI
		Eigen::Vector3d accel_eci = R_ecef2eci*Rsp2cart*accel + accel_central;

		// std::cout << "\n" << "accel ecef (no central): \n" << Rsp2cart*accel;
		// std::cout << "\n" << "accel eci (with central: \n" << accel_eci;

		return accel_eci;

	} //GetGravAccel

} //namespace Util