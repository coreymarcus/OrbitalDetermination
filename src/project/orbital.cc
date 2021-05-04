#include<iostream>
#include<cmath>
#include"../../../../eigen/Eigen/Eigen"
#include<boost/array.hpp>
#include<boost/numeric/odeint.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include "orbital.h"


namespace orbital{

      
      SpaceCraft::SpaceCraft(int d){
        dimensions = d;
        cart.STM = Eigen::MatrixXf::Identity(d,d);
      } // constructor

      Eigen::VectorXf SpaceCraft::getKepl(){}

      Eigen::VectorXf SpaceCraft::getCartInert(){}

      Eigen::VectorXf SpaceCraft::getCartRot(){}
 
      Eigen::MatrixXf LoadDatFile(std::string file, int rows, int cols){

        Eigen::MatrixXf dat = Eigen::MatrixXf::Zero(rows, cols);
        std::ifstream myFile(file);
        if(!myFile.is_open()) throw std::runtime_error("Could not open file");
        std::string line;
	std::string word;
	std::string::size_type sz;
        for (int rowidx = 0; rowidx < rows; ++rowidx) {
          std::getline(myFile, line);
          std::stringstream ss(line);
          for (int colidx = 0; colidx < cols; ++colidx) {
            std::getline(ss, word, ',');
            dat(rowidx,colidx) = std::stod(word, &sz);
          }
        }

        return dat;
      }

      void SpaceCraft::setKepl(Eigen::VectorXf Kepl){
        oe.sma = Kepl[0];
        oe.ecc = Kepl[1];
        oe.incl = Kepl[2];
        oe.RAAN = Kepl[3];
        oe.ome = Kepl[4];
        oe.theta = Kepl[5];
        Kepl2iCart_();
      }

      void SpaceCraft::setCartInert(Eigen::Vector3f ipos, Eigen::Vector3f ivel){
        cart.ipos[0] = ipos[0];
        cart.ipos[1] = ipos[1];
        cart.ipos[2] = ipos[2];
        cart.ivel[0] = ivel[0];
        cart.ivel[1] = ivel[1];
        cart.ivel[2] = ivel[2];
        iCart2Kepl_();
      }

      void SpaceCraft::setCartInertCov(Eigen::MatrixXf covariance){
        if (covariance.size() != dimensions*dimensions) {
          std::cout << "Wrong dimension of the Covariance, it should be " << dimensions << std::endl;
          exit(0);
        }
        cart.Cov = covariance;

      }

      Eigen::MatrixXf SpaceCraft::getQmatrix(double deltaT) {

        Eigen::MatrixXf Q = Eigen::MatrixXf::Zero(dimensions, dimensions);
        Eigen::MatrixXf Gamma = Eigen::MatrixXf::Zero(6,3);
        Gamma.block(0,0,3,3) = Eigen::Matrix3f::Identity(3,3)*(deltaT)*(deltaT)/2.0;
        Gamma.block(3,0,3,3) = Eigen::Matrix3f::Identity(3,3)*std::abs(deltaT);
        Q.block(0,0,6,6) = sigma_Q*sigma_Q*Gamma*Gamma.transpose();
        Q(6,6) = (0.1/5e3)*(0.1/5e3)*std::pow(std::abs(deltaT),2);
        return Q;
      }


      void SpaceCraft::KalmanUpdate(std::vector<double> measurement, int station){
        // Only provide update for radar measurement; only range and range-rate
        // last measurement should be range-rate; before last should be range
        // Expected Measurement take lighttime into account, but update occurs at timestamp of measurement
        std::vector<double> expmeas;
        std::vector<int> indices;
        indices.push_back(0);
        indices.push_back(1);
       
        double lt_corr =  getExpectedMeasurement(station, indices, t, &expmeas); 
        propagate(-lt_corr);


        Eigen::MatrixXf H = getHmatrix(station, t);

        Eigen::Vector2f residuals;
        Eigen::MatrixXf residualsCov;
        residuals(0) = measurement[0] - expmeas[0];
        residuals(1) = measurement[1] - expmeas[1];

        Eigen::Matrix2f R = Eigen::Matrix2f::Zero(2,2);
        R(0,0) = 0.01*0.01;
        R(1,1) = 1e-6*1e-6/4;
        if (station==1) {
          R(0,0) /= 4;
          R(1,1) *= 4;
        }
        //R(0,0) *= 1e30;
        
        bool all_stations = false;
        bool station2 = false;
        bool station1 = true;
        sigma_Q = 100e-9/1e3;
        //sigma_Q = 10e-8;  // station 1
        if ((t > 5e3) ){//and (t < 81220)) {
          sigma_Q = 100e-9/1e3;
        //  sigma_Q = 2e-9;  // station 1
        }
        //if ((t > 15e3) )  // station 1
        //  sigma_Q = 1e-9;
        //}
        Eigen::MatrixXf Q = getQmatrix(t - last_update_time); 
        Eigen::MatrixXf predCov = cart.STM*cart.Cov*cart.STM.transpose() + Q;
        residualsCov = H*predCov*H.transpose() + R;
        //std::cout << t << " , " << std::sqrt(residualsCov(0,0)) <<  " , " << std::sqrt(residualsCov(1,1)) << std::endl;

        Eigen::MatrixXf gain = predCov*H.transpose()*residualsCov.inverse();

        Eigen::VectorXf predState(dimensions);
        predState(0) = cart.ipos(0);
        predState(1) = cart.ipos(1);
        predState(2) = cart.ipos(2);
        predState(3) = cart.ivel(0);
        predState(4) = cart.ivel(1);
        predState(5) = cart.ivel(2);
        predState(6) = env.Cd;
        Eigen::VectorXf postState = predState + gain*residuals;

        cart.ipos(0) = postState(0);
        cart.ipos(1) = postState(1);
        cart.ipos(2) = postState(2);
        cart.ivel(0) = postState(3);
        cart.ivel(1) = postState(4);
        cart.ivel(2) = postState(5);
        env.Cd = postState(6);
        if (env.Cd < 0.5) {
          env.Cd = 0.5;
        } else if (env.Cd > 3.5) {
          env.Cd = 3.5;
        }

        Eigen::MatrixXf I = Eigen::MatrixXf::Identity(dimensions,dimensions);
        cart.Cov = (I - gain*H)*predCov*((I - gain*H)).transpose() + gain*R*gain.transpose();

        cart.STM = Eigen::MatrixXf::Identity(dimensions,dimensions);
        last_update_time = t;
        propagate(lt_corr);

      }
 
      void SpaceCraft::setCartRot(Eigen::Vector3f ipos, Eigen::Vector3f ivel){}

      void SpaceCraft::iCart2Kepl_(){ // convert from Cartesian inertial state to Keplerian state

        // generate versors
        Eigen::Vector3f ii(1.,0.,0.);
        Eigen::Vector3f jj(0.,1.,0.);
        Eigen::Vector3f kk(0.,0.,1.);

        // radius and velocity
        double r = cart.ipos.norm();
        double V = cart.ivel.norm();

        // energy_pm and sma (and period)
        energy_pm = -env.mu/r + 0.5*V*V;
        double R_e_r = env.R_e/r;
        energy_J2 = energy_pm + env.mu/r*env.J_2*R_e_r*R_e_r*(1.5*cart.ipos[2]*cart.ipos[2]/r/r - 0.5);
        oe.sma = -env.mu/2/energy_pm;
        oe.n = 1.0/std::sqrt(oe.sma*oe.sma*oe.sma/env.mu); 
        period = 2.0*M_PI/oe.n;

        // angular momentum
        hh = cart.ipos.cross(cart.ivel);
        h = hh.norm();

        // Eccentricity
        oe.ee = 1./env.mu * (cart.ivel.cross(hh) - env.mu/r*cart.ipos);
        oe.ecc = oe.ee.norm();

        // Inclination
        oe.incl = std::acos(hh[2]/h);

	// Ascendant node vector
	Eigen::Vector3f nn = kk.cross(hh);
        double n = nn.norm();

        // RAAN
        oe.RAAN = std::atan2(nn[1]/n, nn[0]/n);
        if (oe.RAAN < 0.) {
          oe.RAAN += 2*M_PI;
        }

        // Arg of periapsis
        double omeCandidate = std::acos(nn.dot(oe.ee)/oe.ecc/n);
        if (oe.ee[2] >= 0) {
          oe.ome = omeCandidate;
        } else { oe.ome = 2*M_PI - omeCandidate;
        }

        double cosThetaCandidate = (cart.ipos.dot(oe.ee)/r/oe.ecc);
        double thetaCandidate;
        // some if statements to solve potential numeric issues
        if (cosThetaCandidate > 1 - 1e-14) {
          thetaCandidate = 0.;
        } else if ( cosThetaCandidate < - 1 +1e-14) {
          thetaCandidate = M_PI;
        } else { 
          thetaCandidate = std::acos(cosThetaCandidate);
        }

        if (cart.ipos.dot(cart.ivel) >=0) {
          oe.theta = thetaCandidate;
        } else { oe.theta = 2*M_PI - thetaCandidate;
        }
        oe.E = 2.0*std::atan(std::tan(oe.theta/2.0) * std::sqrt((1.0 - oe.ecc)/(1.0 + oe.ecc)));
        if (oe.E<0) {
          oe.E = oe.E + 2.0 * M_PI;
        }
        
        oe.M = oe.E - oe.ecc * std::sin(oe.E);
        
        oe.T_p = t - oe.M/oe.n;

      }

      Eigen::Vector3f SpaceCraft::SunPosition(double JD_UT1, double JD_TDB) { // provides position of the Sun in ECI inertial coordinates
        // Using as input JD_UT1 instead of JD_TDB is not an issue (the Earth is slow moving)
        //
        double T_UT1 = (JD_UT1 - 2451545.)/36525;
        double T_TDB = (JD_TDB - 2451545.)/36525;
        double deg2rad = M_PI/180.;
    
        double lambda_M_sun = (280.460 + 36000.771 * T_UT1)*deg2rad;
        double M_sun = (357.5291092 + 35999.05034*T_TDB)*deg2rad;
        double lambda_ecl = lambda_M_sun + (1.914666471*std::sin(M_sun) 
          + 0.000139589*std::cos(2*M_sun))*deg2rad; 
        double dist_sun = env.AU * (1.000140612 - 0.016708617*std::cos(M_sun)
          - 0.000139589*std::cos(2*M_sun));
        double eps = (23.439291 - 0.0130042*T_TDB)*deg2rad;
        Eigen::Vector3f dirSun;
        dirSun[0] = std::cos(lambda_ecl);
        dirSun[1] = std::cos(eps) * std::sin(lambda_ecl);
        dirSun[2] = std::sin(eps) * std::sin(lambda_ecl);
        return dist_sun * dirSun;

      }

      Eigen::Vector3f SpaceCraft::MoonPosition(double JD_TDB) { // provides position of the Moon in inertial ECI coordinates
        // Using as input JD_UT1 instead of JD_TDB is not an issue (the Moon is slow moving)
        //
        double T_TDB = (JD_TDB - (2451545. - 2400000.5))/36525;
        T_TDB = (JD_TDB - 2451545.)/36525;
        double deg2rad = M_PI/180.;

        double lambda_ecl = (218.32 + 481267.8813*T_TDB + 6.29*std::sin((134.9 + 477198.85*T_TDB)*deg2rad)
          - 1.27*std::sin((259.2 - 413335.38*T_TDB)*deg2rad) + 0.66*std::sin((235.7 + 890534.23*T_TDB)*deg2rad)
          + 0.21*std::sin((269.9 + 954397.70*T_TDB)*deg2rad) - 0.19*std::sin((357.5 + 35999.05*T_TDB)*deg2rad)
          - 0.11*std::sin((186.6 + 966404*T_TDB)*deg2rad)
          )*deg2rad; //radians
        double phi_ecl = (5.13*std::sin((93.3 + 483202.03*T_TDB)*deg2rad)
          + 0.28*std::sin((228.2 + 960400.87*T_TDB)*deg2rad) - 0.28*std::sin((318.3 + 6003.18*T_TDB)*deg2rad)
          - 0.17*std::sin((217.6 - 407332.2*T_TDB)*deg2rad)
          )*deg2rad; //radians
        double cal_P = (0.9508 + 0.0518*std::cos((134.9 + 477198.85*T_TDB)*deg2rad)
          + 0.0095*std::cos((259.2 - 413335.38*T_TDB)*deg2rad) + 0.0078*std::cos((235.7 + 890534.23*T_TDB)*deg2rad)
          + 0.0028*std::cos((269.9 + 954397.70*T_TDB)*deg2rad)
          )*deg2rad; //radians
        double eps_bar = (23.439291 - 0.0130042*T_TDB - 1.64e-7*T_TDB*T_TDB + 5.04e-7*T_TDB*T_TDB*T_TDB
          )*deg2rad; // radians
        
        double distMoon = 1./std::sin(cal_P)*env.R_e; // this is going to be around 380 000 km
        Eigen::Vector3f dirMoon;
        dirMoon[0] = std::cos(phi_ecl) * std::cos(lambda_ecl);
        dirMoon[1] = std::cos(eps_bar) * std::cos(phi_ecl) * std::sin(lambda_ecl) 
                   - std::sin(eps_bar) * std::sin(phi_ecl);
        dirMoon[2] = std::sin(eps_bar) * std::cos(phi_ecl) * std::sin(lambda_ecl) 
                   + std::cos(eps_bar) * std::sin(phi_ecl);
        return distMoon * dirMoon;

      }
      Eigen::Vector3f SpaceCraft::Cart2spher_vec_(Eigen::Vector3f cartvec){
        // convert from Cartesian coord to spherical coord
        // sphervec(0) is r
        // sphervec(1) is longitute (-pi, pi)
        // sphervec(2) is latitude (-pi/2,pi/2)
        Eigen::Vector3f sphervec;
        sphervec(0) = cartvec.norm();
        sphervec(1) = std::atan2(cartvec(1), cartvec(0));
        sphervec(2) = std::asin(cartvec(2)/sphervec(0));
        return sphervec;
      }  

      void SpaceCraft::Kepl2iCart_(){ // convert from Keplerian state to Cartesian inertial state

        // Calculating parameter p and radius r
        double p = oe.sma*(1-oe.ecc*oe.ecc);
        double r = p/(1+oe.ecc*std::cos(oe.theta));

        // Calculating xi and eta
        double xi = r*std::cos(oe.theta);
        double eta = r*std::sin(oe.theta);
        Eigen::VectorXf xieta(2);
        xieta << xi, eta;

        // Caclulating l, m and n
        double l1 = std::cos(oe.RAAN)*std::cos(oe.ome) - std::sin(oe.RAAN)*std::sin(oe.ome)*std::cos(oe.incl);
        double l2 = -std::cos(oe.RAAN)*std::sin(oe.ome) - std::sin(oe.RAAN)*std::cos(oe.ome)*std::cos(oe.incl);
        double m1 = std::sin(oe.RAAN)*std::cos(oe.ome) + std::cos(oe.RAAN)*std::sin(oe.ome)*std::cos(oe.incl);
        double m2 = - std::sin(oe.RAAN)*std::sin(oe.ome) + std::cos(oe.RAAN)*std::cos(oe.ome)*std::cos(oe.incl);
        double n1 = std::sin(oe.ome)*std::sin(oe.incl);
        double n2 = std::cos(oe.ome)*std::sin(oe.incl);

        // Building M matrix
        Eigen::MatrixXf M(3,2);
        M << l1, l2,
             m1, m2,
             n1, n2;
        Eigen::VectorXf rmx = M*xieta;
        cart.ipos[0] = rmx[0];
        cart.ipos[1] = rmx[1];
        cart.ipos[2] = rmx[2];

        double H = std::sqrt(env.mu*p);
	cart.ivel[0] = env.mu/H*(-l1*std::sin(oe.theta)+l2*(oe.ecc+std::cos(oe.theta)));
        cart.ivel[1] = env.mu/H*(-m1*std::sin(oe.theta)+m2*(oe.ecc+std::cos(oe.theta)));
        cart.ivel[2] = env.mu/H*(-n1*std::sin(oe.theta)+n2*(oe.ecc+std::cos(oe.theta)));
        
      }

      void SpaceCraft::setSTMprop() {
        STM_PROP_ = true;
      }





      state_type SpaceCraft::STM_dot2 (const state_type X , const double t ){
        // Provides derivative of STM with J2 and J3
        // Input (and output) are (d*d+1) vectors

        state_type dXdt(dimensions*(dimensions+1),0.0);

        double x = X[0];        
        double y = X[1];        
        double z = X[2];        
        double Vx = X[3];        
        double Vy = X[4];        
        double Vz = X[5];        

        int d = dimensions;
 
        double R_e = env.R_e;
        double mu = env.mu;

        double H;
        double omega;
        double ball;
        double rho_0;
        double h_0;
        double Cd;


        if (EXP_DRAGCANNON_) {
          H = env.H;
          omega = env.omega;
          ball = env.mass/env.Cd/env.Area;
          rho_0 = env.rho_0;
          h_0 = env.r_0 - R_e;
          Cd = env.Cd;
        } else { ball = 100;
          omega = 0.;
          rho_0 = 0;
          H = 10.;
          h_0 = 100;
          ball = 100.;
          Cd = 2.0;
        }
        double J2;
        double J3;
        if (LEG_GRAV_FIELD_) {
          J2 = - env.Cgrav(2,0);
          J3 = - env.Cgrav(3,0);
        } else { J2 = env.J_2;
          J3 = env.J_3;
        }
        
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(d,d);
        Eigen::MatrixXf Phi = Eigen::MatrixXf::Zero(d,d);

        double t1 = x * x;
        double t2 = y * y;
        double t3 = z * z;
        double t4 = t1 + t2 + t3;
        double t5 = t4 * t4;
        double t6 = t5 * t5;
        double t8 = sqrt(t4);
        double t10 = 0.1e1 / t8 / t4 / t6;
        double t11 = H * ball;
        double t12 = t1 * t1;
        double t13 = t12 * t12;
        double t15 = t13 * mu * t11;
        double t17 = H * J2;
        double t18 = R_e * R_e;
        double t19 = ball * t18;
        double t21 = mu * t19 * t17;
        double t23 = t2 * mu;
        double t24 = t23 * t11;
        double t27 = t3 * mu * t11;
        double t30 = t1 * t12;
        double t32 = t2 * t2;
        double t33 = t32 * mu;
        double t34 = t33 * t11;
        double t40 = t18 * t17;
        double t41 = mu * ball;
        double t42 = t3 * t41;
        double t43 = t42 * t40;
        double t45 = J3 * H;
        double t46 = R_e * t18;
        double t47 = t46 * t45;
        double t48 = z * t41;
        double t49 = t48 * t47;
        double t51 = t3 * t3;
        double t53 = t51 * mu * t11;
        double t54 = 0.6666666666e-1 * t53;
        double t57 = t2 * t32;
        double t58 = t57 * mu;
        double t59 = t58 * t11;
        double t69 = t51 * t41;
        double t70 = t69 * t40;
        double t72 = z * t3;
        double t73 = t72 * t41;
        double t74 = t73 * t47;
        double t76 = t3 * t51;
        double t78 = t76 * mu * t11;
        double t82 = t32 * t32;
        double t84 = t82 * mu * t11;
        double t100 = t51 * t51;
        double t102 = t100 * mu * t11;
        double t105 = t76 * t41 * t40;
        double t107 = z * t51;
        double t108 = t107 * t41;
        double t109 = t108 * t47;
        double t111 = 0.4444444444e-1 * t15 + t30 * (0.1333333333e0 * t21 + 0.1111111111e0 * t24 + 0.1111111111e0 * t27) + t12 * (0.6666666666e-1 * t34 + t2 * (0.1333333333e0 * t27 + 0.2333333333e0 * t21) - 0.7666666666e0 * t43 + 0.9999999999e0 * t49 + t54) + t1 * (-0.2222222222e-1 * t59 + t32 * (-0.6666666666e-1 * t27 + 0.6666666666e-1 * t21) + t2 * (-t54 - 0.6999999999e0 * t43 + 0.8333333332e0 * t49) - 0.7666666666e0 * t70 - 0.2277777778e1 * t74 - 0.2222222222e-1 * t78) - 0.2222222222e-1 * t84 + t57 * (-0.3333333333e-1 * t21 - 0.8888888888e-1 * t27) + t32 * (0.6666666666e-1 * t43 - 0.1666666666e0 * t49 - 0.1333333333e0 * t53) + t2 * (0.2333333333e0 * t70 + 0.5555555555e-1 * t74 - 0.8888888888e-1 * t78) - 0.2222222222e-1 * t102 + 0.1333333333e0 * t105 + 0.2222222222e0 * t109;
        double t113 = omega * omega;
        double t114 = t113 * (t1 + t2);
        double t115 = y * Vx;
        double t117 = x * Vy;
        double t121 = Vx * Vx;
        double t122 = Vy * Vy;
        double t123 = Vz * Vz;
        double t125 = sqrt(t114 + omega * (0.2e1 * t115 - 0.2e1 * t117) + t121 + t122 + t123);
        double t127 = Vx * H;
        double t130 = omega * t113;
        double t131 = t130 * H;
        double t136 = x * t1;
        double t137 = t136 * t13;
        double t139 = H * Vy;
        double t140 = y * t113;
        double t143 = Vy * omega;
        double t144 = t143 * t127;
        double t148 = t1 * t13;
        double t150 = y * t2;
        double t153 = t3 * t113;
        double t154 = t153 * t127;
        double t156 = t2 * t113;
        double t159 = t3 * y;
        double t164 = x * t13;
        double t166 = t3 * t143;
        double t167 = t166 * t127;
        double t169 = t3 * t140;
        double t172 = t2 * t143;
        double t175 = t150 * t113;
        double t181 = t51 * y;
        double t184 = t3 * t150;
        double t187 = y * t32;
        double t190 = t3 * t156;
        double t193 = t32 * t113;
        double t196 = t51 * t113;
        double t197 = t196 * t127;
        double t201 = t136 * t12;
        double t203 = t187 * t113;
        double t206 = t32 * t143;
        double t209 = t51 * t140;
        double t212 = t3 * t175;
        double t215 = t51 * t143;
        double t216 = t215 * t127;
        double t218 = Vy * t127;
        double t219 = t2 * omega;
        double t220 = t3 * t219;
        double t226 = t3 * t193;
        double t229 = t51 * t156;
        double t232 = t76 * t113;
        double t233 = t232 * t127;
        double t235 = t51 * t150;
        double t238 = t150 * t32;
        double t241 = t3 * t187;
        double t244 = t76 * y;
        double t247 = t57 * t113;
        double t252 = x * t12;
        double t254 = t76 * t143;
        double t255 = t254 * t127;
        double t257 = t57 * t143;
        double t260 = t76 * t140;
        double t263 = t51 * t175;
        double t266 = t3 * t203;
        double t269 = t238 * t113;
        double t272 = t51 * t219;
        double t275 = t32 * omega;
        double t276 = t3 * t275;
        double t282 = t100 * y;
        double t285 = t3 * t247;
        double t288 = t51 * t193;
        double t291 = y * t82;
        double t294 = t76 * t156;
        double t297 = t100 * t113;
        double t298 = t297 * t127;
        double t300 = t51 * t187;
        double t303 = t82 * t113;
        double t306 = t76 * t150;
        double t309 = t3 * t238;
        double t315 = t82 * t143;
        double t318 = t100 * t143;
        double t319 = t318 * t127;
        double t321 = t3 * t269;
        double t324 = t51 * t203;
        double t327 = t76 * t175;
        double t330 = t100 * t140;
        double t333 = t291 * t113;
        double t336 = t51 * t275;
        double t339 = t76 * t219;
        double t342 = t57 * omega;
        double t343 = t3 * t342;
        double t349 = t100 * t150;
        double t352 = t51 * t247;
        double t355 = t100 * t156;
        double t358 = t3 * t303;
        double t361 = t76 * t193;
        double t364 = t150 * t82;
        double t367 = t51 * t238;
        double t370 = t76 * t187;
        double t373 = t3 * t100;
        double t374 = t373 * y;
        double t377 = t2 * t82;
        double t378 = t377 * t113;
        double t381 = t3 * t291;
        double t384 = t373 * t113;
        double t385 = t384 * t127;
        double t387 = -0.5555555555e-1 * t349 * t131 - 0.1111111111e0 * t352 * t127 - 0.5555555555e-1 * t355 * t127 - 0.5555555555e-1 * t358 * t127 - 0.1111111111e0 * t361 * t127 - 0.1111111111e-1 * t364 * t131 - 0.1111111111e0 * t367 * t131 - 0.1111111111e0 * t370 * t131 - 0.1111111111e-1 * t374 * t131 - 0.1111111111e-1 * t378 * t127 - 0.5555555555e-1 * t381 * t131 - 0.1111111111e-1 * t385;
        double t390 = t373 * t140;
        double t393 = t373 * t143;
        double t394 = t393 * t127;
        double t396 = t76 * t275;
        double t399 = t51 * t342;
        double t402 = t100 * t219;
        double t405 = t82 * omega;
        double t406 = t3 * t405;
        double t409 = t364 * t113;
        double t412 = t51 * t269;
        double t415 = t377 * t143;
        double t418 = t76 * t203;
        double t421 = t100 * t175;
        double t424 = t3 * t333;
        double t427 = 0.1111111111e-1 * t390 * t139 + 0.1111111111e-1 * t394 + 0.1111111111e0 * t396 * t218 + 0.1111111111e0 * t399 * t218 + 0.5555555555e-1 * t402 * t218 + 0.5555555555e-1 * t406 * t218 + 0.1111111111e-1 * t409 * t139 + 0.1111111111e0 * t412 * t139 + 0.1111111111e-1 * t415 * t127 + 0.1111111111e0 * t418 * t139 + 0.5555555555e-1 * t421 * t139 + 0.5555555555e-1 * t424 * t139;
        double t429 = t137 * rho_0 * (-0.1111111111e-1 * t113 * t127 - 0.1111111111e-1 * y * t131) + t148 * rho_0 * (0.1111111111e-1 * t140 * t139 + 0.1111111111e-1 * t144) + t164 * rho_0 * (-0.5555555555e-1 * t150 * t131 - 0.5555555555e-1 * t154 - 0.5555555555e-1 * t156 * t127 - 0.5555555555e-1 * t159 * t131) + t13 * rho_0 * (0.5555555555e-1 * t167 + 0.5555555555e-1 * t169 * t139 + 0.5555555555e-1 * t172 * t127 + 0.5555555555e-1 * t175 * t139) + t201 * rho_0 * (-0.1111111111e0 * t181 * t131 - 0.2222222222e0 * t184 * t131 - 0.1111111111e0 * t187 * t131 - 0.2222222222e0 * t190 * t127 - 0.1111111111e0 * t193 * t127 - 0.1111111111e0 * t197) + t30 * rho_0 * (0.1111111111e0 * t203 * t139 + 0.1111111111e0 * t206 * t127 + 0.1111111111e0 * t209 * t139 + 0.2222222222e0 * t212 * t139 + 0.1111111111e0 * t216 + 0.2222222222e0 * t220 * t218) + t252 * rho_0 * (-0.3333333333e0 * t226 * t127 - 0.3333333333e0 * t229 * t127 - 0.1111111111e0 * t233 - 0.3333333333e0 * t235 * t131 - 0.1111111111e0 * t238 * t131 - 0.3333333333e0 * t241 * t131 - 0.1111111111e0 * t244 * t131 - 0.1111111111e0 * t247 * t127) + t12 * rho_0 * (0.1111111111e0 * t255 + 0.1111111111e0 * t257 * t127 + 0.1111111111e0 * t260 * t139 + 0.3333333333e0 * t263 * t139 + 0.3333333333e0 * t266 * t139 + 0.1111111111e0 * t269 * t139 + 0.3333333333e0 * t272 * t218 + 0.3333333333e0 * t276 * t218) + t136 * rho_0 * (-0.5555555555e-1 * t282 * t131 - 0.2222222222e0 * t285 * t127 - 0.3333333333e0 * t288 * t127 - 0.5555555555e-1 * t291 * t131 - 0.2222222222e0 * t294 * t127 - 0.5555555555e-1 * t298 - 0.3333333333e0 * t300 * t131 - 0.5555555555e-1 * t303 * t127 - 0.2222222222e0 * t306 * t131 - 0.2222222222e0 * t309 * t131) + t1 * rho_0 * (0.5555555555e-1 * t315 * t127 + 0.5555555555e-1 * t319 + 0.2222222222e0 * t321 * t139 + 0.3333333333e0 * t324 * t139 + 0.2222222222e0 * t327 * t139 + 0.5555555555e-1 * t330 * t139 + 0.5555555555e-1 * t333 * t139 + 0.3333333333e0 * t336 * t218 + 0.2222222222e0 * t339 * t218 + 0.2222222222e0 * t343 * t218) + x * rho_0 * t387 + rho_0 * t427;
        double t431 = t113 * Vx;
        double t437 = t252 * t13;
        double t439 = Vy * Vx;
        double t440 = omega * t439;
        double t442 = t113 * Vy;
        double t447 = t12 * t13;
        double t453 = t3 * t130;
        double t456 = 0.1111111111e-1 * t122;
        double t457 = 0.1111111111e-1 * t123;
        double t458 = 0.3333333333e-1 * t121 + t456 + t457;
        double t460 = 0.5555555555e-1 * t453 + omega * t458;
        double t462 = t3 * t431;
        double t464 = Vx * t121;
        double t465 = 0.1111111111e-1 * t464;
        double t467 = Vx * (t456 + t457);
        double t473 = t3 * omega;
        double t474 = t473 * t439;
        double t488 = 0.5555555555e-1 * t122;
        double t489 = 0.5555555555e-1 * t123;
        double t491 = t488 + t489 + 0.1666666666e0 * t121;
        double t493 = 0.2777777778e0 * t453 + omega * t491;
        double t496 = 0.5555555555e-1 * t464;
        double t498 = Vx * (t488 + t489);
        double t501 = t51 * t130;
        double t505 = 0.1111111111e0 * t501 + omega * t3 * t491;
        double t507 = t51 * t431;
        double t509 = t496 + t498;
        double t510 = t3 * t509;
        double t524 = t51 * omega;
        double t525 = t524 * t439;
        double t535 = 0.1111111111e0 * t123;
        double t537 = 0.1111111111e0 * t122;
        double t538 = t535 + 0.3333333333e0 * t121 + t537;
        double t540 = 0.5555555555e0 * t453 + omega * t538;
        double t543 = 0.1111111111e0 * t464;
        double t545 = Vx * (t535 + t537);
        double t550 = 0.2222222222e0 * t122;
        double t551 = 0.2222222222e0 * t123;
        double t552 = 0.6666666666e0 * t121 + t550 + t551;
        double t555 = 0.4444444444e0 * t501 + omega * t3 * t552;
        double t561 = 0.2222222222e0 * t464 + Vx * (t550 + t551);
        double t562 = t3 * t561;
        double t565 = t76 * t130;
        double t569 = 0.1111111111e0 * t565 + omega * t51 * t538;
        double t571 = t76 * t431;
        double t573 = t543 + t545;
        double t574 = t51 * t573;
        double t586 = t76 * omega;
        double t587 = t586 * t439;
        double t608 = 0.3333333333e0 * t122;
        double t609 = 0.3333333333e0 * t123;
        double t610 = 0.9999999999e0 * t121 + t608 + t609;
        double t619 = 0.3333333333e0 * t464 + Vx * (t608 + t609);
        double t626 = 0.3333333333e0 * t565 + omega * t51 * t610;
        double t629 = t51 * t619;
        double t632 = t100 * t130;
        double t636 = 0.5555555555e-1 * t632 + omega * t76 * t538;
        double t638 = t100 * t431;
        double t640 = t76 * t573;
        double t641 = 0.1666666666e0 * t291 * t130 + 0.3888888888e0 * t82 * t431 + t238 * t540 + t57 * (0.1222222222e1 * t462 + t543 + t545) + t187 * (0.6666666666e0 * t501 + omega * t3 * t610) + t32 * (0.1333333333e1 * t507 + t3 * t619) + t150 * t626 + t2 * (0.5555555555e0 * t571 + t629) + y * t636 + 0.5555555555e-1 * t638 + t640;
        double t656 = t100 * omega;
        double t657 = t656 * t439;
        double t693 = t373 * t130;
        double t697 = 0.1111111111e-1 * t693 + omega * t100 * t491;
        double t699 = t373 * t431;
        double t701 = t100 * t509;
        double t702 = 0.6666666666e-1 * t364 * t130 + 0.1777777778e0 * t377 * t431 + t291 * t493 + t82 * (0.7222222222e0 * t462 + t496 + t498) + t238 * t555 + t57 * (0.1111111111e1 * t507 + t562) + t187 * t626 + t32 * (0.7777777777e0 * t571 + t629) + t150 * (0.1111111111e0 * t632 + omega * t76 * t552) + t2 * (0.2222222222e0 * t638 + t76 * t561) + y * t697 + 0.1111111111e-1 * t699 + t701;
        double t717 = t373 * omega;
        double t718 = t717 * t439;
        double t722 = t377 * omega;
        double t731 = -0.1111111111e0 * t406 * t439 - 0.2222222222e0 * t399 * t439 - 0.2222222222e0 * t396 * t439 - 0.1111111111e0 * t402 * t439 - 0.2222222222e0 * t367 * t442 - 0.2222222222e-1 * t364 * t442 - 0.2222222222e-1 * t718 - 0.1111111111e0 * t381 * t442 - 0.2222222222e-1 * t722 * t439 - 0.2222222222e0 * t370 * t442 - 0.1111111111e0 * t349 * t442 - 0.2222222222e-1 * t374 * t442;
        double t734 = t187 * t82;
        double t737 = t32 * t82;
        double t761 = omega * y;
        double t765 = 0.1111111111e-1 * t734 * t130 + 0.3333333333e-1 * t737 * t431 + t364 * t460 + t377 * (0.1666666666e0 * t462 + t465 + t467) + t291 * t505 + t82 * (0.3333333333e0 * t507 + t510) + t238 * t569 + t57 * (0.3333333333e0 * t571 + t574) + t187 * t636 + t32 * (0.1666666666e0 * t638 + t640) + t150 * t697 + t2 * (0.3333333333e-1 * t699 + t701) + t761 * t373 * t458 + t373 * (t465 + t467);
        double t768 = t8 * t429 + t437 * rho_0 * (0.1111111111e-1 * t431 + 0.1111111111e-1 * y * t130) + t447 * rho_0 * (-0.2222222222e-1 * t440 - 0.2222222222e-1 * y * t442) + t137 * rho_0 * (0.6666666666e-1 * t150 * t130 + 0.8888888888e-1 * t2 * t431 + y * t460 + 0.5555555555e-1 * t462 + t465 + t467) + t148 * rho_0 * (-0.1111111111e0 * t219 * t439 - 0.1111111111e0 * t474 - 0.1111111111e0 * t159 * t442 - 0.1111111111e0 * t150 * t442) + t164 * rho_0 * (0.1666666666e0 * t187 * t130 + 0.2777777778e0 * t32 * t431 + t150 * t493 + t2 * (0.3888888888e0 * t462 + t496 + t498) + y * t505 + 0.1111111111e0 * t507 + t510) + t13 * rho_0 * (-0.4444444444e0 * t220 * t439 - 0.2222222222e0 * t187 * t442 - 0.4444444444e0 * t184 * t442 - 0.2222222222e0 * t181 * t442 - 0.2222222222e0 * t275 * t439 - 0.2222222222e0 * t525) + t201 * rho_0 * (0.2222222222e0 * t238 * t130 + 0.4444444444e0 * t57 * t431 + t187 * t540 + t32 * (0.9999999999e0 * t462 + t543 + t545) + t150 * t555 + t2 * (0.6666666666e0 * t507 + t562) + y * t569 + 0.1111111111e0 * t571 + t574) + t30 * rho_0 * (-0.6666666666e0 * t276 * t439 - 0.6666666666e0 * t241 * t442 - 0.2222222222e0 * t244 * t442 - 0.2222222222e0 * t342 * t439 - 0.2222222222e0 * t587 - 0.2222222222e0 * t238 * t442 - 0.6666666666e0 * t235 * t442 - 0.6666666666e0 * t272 * t439) + t252 * rho_0 * t641 + t12 * rho_0 * (-0.6666666666e0 * t336 * t439 - 0.4444444444e0 * t339 * t439 - 0.6666666666e0 * t300 * t442 - 0.4444444444e0 * t309 * t442 - 0.1111111111e0 * t282 * t442 - 0.1111111111e0 * t405 * t439 - 0.1111111111e0 * t657 - 0.4444444444e0 * t306 * t442 - 0.4444444444e0 * t343 * t439 - 0.1111111111e0 * t291 * t442) + t136 * rho_0 * t702 + t1 * rho_0 * t731 + x * rho_0 * t765;
        double t770 = 0.1e1 / H;
        double t772 = exp(t770 * (-t8 + R_e + h_0));
        double t779 = sqrt(t114 + 0.2e1 * omega * (t115 - t117) + t121 + t122 + t123);
        double t780 = 0.1e1 / t779;
        double t782 = 0.1e1 / ball;
        double t783 = t782 * t770 * t780;
        double t786 = x * mu;
        double t791 = t136 * mu * t11;
        double t795 = 0.1428571429e0 * t21 + 0.1714285714e0 * t27;
        double t800 = t252 * mu * t11;
        double t804 = 0.2857142858e0 * t21 + 0.3428571429e0 * t27;
        double t809 = 0.1714285714e0 * t53 - 0.7142857144e0 * t43 + 0.1000000000e1 * t49;
        double t813 = t201 * mu;
        double t814 = t813 * t11;
        double t821 = 0.5714285715e-1 * t78 - 0.8571428572e0 * t70 - 0.2000000000e1 * t74;
        double t827 = t136 * t113;
        double t828 = t100 * t827;
        double t831 = t164 * t113;
        double t832 = t3 * t831;
        double t835 = t252 * t113;
        double t836 = t76 * t835;
        double t839 = t201 * t113;
        double t840 = t51 * t839;
        double t843 = x * t113;
        double t844 = t373 * t843;
        double t847 = t76 * t131;
        double t848 = 0.9523809525e-1 * t847;
        double t849 = t121 * H;
        double t851 = 0.9523809525e-1 * t122;
        double t852 = 0.9523809525e-1 * t123;
        double t853 = -t851 - t852;
        double t855 = -0.1904761905e0 * t849 + H * t853;
        double t857 = omega * t51 * t855;
        double t864 = t3 * t131;
        double t865 = 0.2857142858e0 * t864;
        double t867 = 0.4761904762e-1 * t122;
        double t868 = 0.4761904762e-1 * t123;
        double t869 = -t867 - t868;
        double t871 = -0.9523809525e-1 * t849 + H * t869;
        double t872 = omega * t871;
        double t875 = t3 * t839;
        double t878 = t51 * t131;
        double t879 = 0.4761904762e0 * t878;
        double t881 = 0.1904761905e0 * t123;
        double t882 = 0.1904761905e0 * t122;
        double t883 = -t881 - t882;
        double t885 = -0.3809523810e0 * t849 + H * t883;
        double t887 = omega * t3 * t885;
        double t890 = t51 * t835;
        double t893 = 0.3809523810e0 * t847;
        double t895 = 0.2857142858e0 * t123;
        double t896 = 0.2857142858e0 * t122;
        double t897 = -t895 - t896;
        double t899 = -0.5714285715e0 * t849 + H * t897;
        double t901 = omega * t51 * t899;
        double t904 = t76 * t827;
        double t907 = t100 * t131;
        double t908 = 0.1428571429e0 * t907;
        double t913 = t100 * t843;
        double t916 = t373 * t131;
        double t917 = 0.1904761905e-1 * t916;
        double t919 = omega * t100 * t871;
        double t920 = -0.6666666668e-1 * t148 * t131 + 0.9523809525e-1 * t831 * t139 + t13 * (-t865 + t872) + 0.3809523810e0 * t875 * t139 + t30 * (-t879 + t887) + 0.5714285715e0 * t890 * t139 + t12 * (-t893 + t901) + 0.3809523810e0 * t904 * t139 + t1 * (omega * t76 * t885 - t908) + 0.9523809525e-1 * t913 * t139 - t917 + t919;
        double t922 = t1 * t113;
        double t923 = t3 * t922;
        double t927 = t12 * t113;
        double t932 = t76 * t927;
        double t935 = t30 * t113;
        double t936 = t51 * t935;
        double t939 = t13 * t113;
        double t940 = t3 * t939;
        double t943 = t100 * t922;
        double t947 = t148 * t113;
        double t956 = 0.4285714286e0 * t864;
        double t959 = t3 * t843;
        double t962 = 0.1904761905e0 * t878;
        double t964 = omega * t3 * t871;
        double t967 = 0.4761904762e-1 * t907;
        double t969 = omega * t76 * t855;
        double t972 = 0.9523809525e-2 * t916;
        double t975 = 0.9523809525e-1 * t828 * t139 + 0.9523809525e-1 * t832 * t139 + 0.1904761905e0 * t836 * t139 + 0.1904761905e0 * t840 * t139 + 0.1904761905e-1 * t844 * t139 + t30 * (-t848 + t857) + t2 * t920 + t238 * (-0.7619047620e0 * t923 * t127 - 0.3809523810e0 * t197 - 0.3809523810e0 * t927 * t127) + y * (-0.3809523810e0 * t932 * t127 - 0.3809523810e0 * t936 * t127 - 0.1904761905e0 * t940 * t127 - 0.1904761905e0 * t943 * t127 - 0.3809523810e-1 * t385 - 0.3809523810e-1 * t947 * t127) + t82 * (-0.2380952381e0 * t12 * t131 + 0.9523809525e-1 * t827 * t139 + t1 * (-t956 + t872) + 0.9523809525e-1 * t959 * t139 - t962 + t964) + t12 * (-t967 + t969) + t1 * (-t972 + t919);
        double t976 = t137 * t113;
        double t977 = t976 * t139;
        double t979 = t409 * t127;
        double t985 = 0.7619047620e0 * t864;
        double t986 = omega * t855;
        double t989 = t3 * t827;
        double t992 = 0.6666666668e0 * t878;
        double t995 = t51 * t843;
        double t998 = 0.1904761905e0 * t847;
        double t1005 = 0.6666666668e0 * t864;
        double t1008 = t3 * t835;
        double t1011 = 0.8571428572e0 * t878;
        double t1016 = t51 * t827;
        double t1019 = 0.4761904762e0 * t847;
        double t1022 = t76 * t843;
        double t1025 = 0.9523809525e-1 * t907;
        double t1028 = t3 * t927;
        double t1031 = t51 * t922;
        double t1043 = 0.9523809525e-1 * t864;
        double t1045 = 0.9523809525e-2 * t122;
        double t1046 = 0.9523809525e-2 * t123;
        double t1047 = -t1045 - t1046;
        double t1049 = -0.1904761905e-1 * t849 + H * t1047;
        double t1050 = omega * t1049;
        double t1058 = t76 * t922;
        double t1061 = t3 * t935;
        double t1064 = t51 * t927;
        double t1072 = 0.4761904762e-1 * t864;
        double t1075 = 0.9523809525e-1 * t878;
        double t1078 = t737 * t131;
        double t1080 = t447 * t131;
        double t1084 = 0.1904761905e-1 * t977 - 0.3809523810e-1 * t979 + t57 * (-0.2857142858e0 * t30 * t131 + 0.1904761905e0 * t835 * t139 + t12 * (-t985 + t986) + 0.3809523810e0 * t989 * t139 + t1 * (-t992 + t887) + 0.1904761905e0 * t995 * t139 - t998 + t857) + t32 * (-0.1904761905e0 * t13 * t131 + 0.1904761905e0 * t839 * t139 + t30 * (-t1005 + t986) + 0.5714285715e0 * t1008 * t139 + t12 * (omega * t3 * t899 - t1011) + 0.5714285715e0 * t1016 * t139 + t1 * (-t1019 + t901) + 0.1904761905e0 * t1022 * t139 - t1025 + t969) + t187 * (-0.1142857143e1 * t1028 * t127 - 0.1142857143e1 * t1031 * t127 - 0.3809523810e0 * t233 - 0.3809523810e0 * t935 * t127) + t377 * (-0.1047619048e0 * t1 * t131 + 0.1904761905e-1 * t843 * t139 - t1043 + t1050) + t291 * (-0.1904761905e0 * t154 - 0.1904761905e0 * t922 * t127) + t150 * (-0.7619047620e0 * t1058 * t127 - 0.7619047620e0 * t1061 * t127 - 0.1142857143e1 * t1064 * t127 - 0.1904761905e0 * t298 - 0.1904761905e0 * t939 * t127) + t148 * (-t1072 + t1050) + t13 * (-t1075 + t964) - 0.1904761905e-1 * t1078 - 0.9523809525e-2 * t1080 + omega * t373 * t1049;
        double t1097 = 0.4761904762e-1 * t453;
        double t1099 = 0.2857142858e-1 * t121 + t1045 + t1046;
        double t1100 = omega * t1099;
        double t1105 = omega * x;
        double t1109 = 0.9523809525e-2 * t464;
        double t1110 = -Vx * t1047;
        double t1117 = 0.2380952381e0 * t453;
        double t1119 = t867 + 0.1428571429e0 * t121 + t868;
        double t1121 = omega * t1119 + t1117;
        double t1123 = t3 * x;
        double t1126 = 0.9523809525e-1 * t501;
        double t1128 = omega * t3 * t1119;
        double t1133 = t136 * omega;
        double t1137 = 0.4761904762e-1 * t464;
        double t1138 = -Vx * t869;
        double t1145 = t1137 + t1138;
        double t1146 = t3 * t1145;
        double t1153 = 0.4761904762e0 * t453;
        double t1154 = 0.2857142858e0 * t121;
        double t1155 = t1154 + t851 + t852;
        double t1157 = omega * t1155 + t1153;
        double t1159 = t3 * t136;
        double t1162 = 0.3809523810e0 * t501;
        double t1164 = t882 + t881 + 0.5714285715e0 * t121;
        double t1167 = omega * t3 * t1164 + t1162;
        double t1169 = t51 * x;
        double t1172 = 0.9523809525e-1 * t565;
        double t1174 = omega * t51 * t1155;
        double t1179 = t252 * omega;
        double t1183 = 0.9523809525e-1 * t464;
        double t1184 = -Vx * t853;
        double t1193 = 0.1904761905e0 * t464 - Vx * t883;
        double t1194 = t3 * t1193;
        double t1201 = t1183 + t1184;
        double t1202 = t51 * t1201;
        double t1210 = t3 * t252;
        double t1213 = 0.5714285715e0 * t501;
        double t1215 = 0.8571428572e0 * t121 + t896 + t895;
        double t1220 = t51 * t136;
        double t1223 = 0.2857142858e0 * t565;
        double t1226 = omega * t51 * t1215 + t1223;
        double t1228 = t76 * x;
        double t1231 = 0.4761904762e-1 * t632;
        double t1233 = omega * t76 * t1155;
        double t1238 = t201 * omega;
        double t1250 = 0.2857142858e0 * t464 - Vx * t897;
        double t1258 = t51 * t1250;
        double t1265 = t76 * t1201;
        double t1273 = t3 * t201;
        double t1277 = t51 * t252;
        double t1281 = t76 * t136;
        double t1284 = 0.9523809525e-1 * t632;
        double t1289 = t100 * x;
        double t1292 = 0.9523809525e-2 * t693;
        double t1294 = omega * t100 * t1119;
        double t1295 = 0.5714285715e-1 * t148 * t130 - 0.9523809525e-1 * t164 * t442 + t13 * t1121 - 0.3809523810e0 * t1273 * t442 + t30 * t1167 - 0.5714285715e0 * t1277 * t442 + t12 * t1226 - 0.3809523810e0 * t1281 * t442 + t1 * (omega * t76 * t1164 + t1284) - 0.9523809525e-1 * t1289 * t442 + t1292 + t1294;
        double t1299 = t164 * omega;
        double t1328 = t100 * t1145;
        double t1329 = 0.7619047620e-1 * t148 * t431 - 0.9523809525e-1 * t1299 * t439 + t13 * (0.3333333334e0 * t462 + t1137 + t1138) - 0.3809523810e0 * t3 * t1238 * t439 + t30 * (0.5714285715e0 * t507 + t1194) - 0.5714285715e0 * t51 * t1179 * t439 + t12 * (0.4761904762e0 * t571 + t1258) - 0.3809523810e0 * t76 * t1133 * t439 + t1 * (0.1904761905e0 * t638 + t76 * t1193) - 0.9523809525e-1 * t100 * t1105 * t439 + 0.2857142858e-1 * t699 + t1328;
        double t1337 = t3 * t164;
        double t1342 = t51 * t201;
        double t1347 = t76 * t252;
        double t1352 = t100 * t136;
        double t1357 = t373 * x;
        double t1362 = 0.9523809525e-2 * t447 * t130 - 0.1904761905e-1 * t137 * t442 + t148 * (t1097 + t1100) - 0.9523809525e-1 * t1337 * t442 + t13 * (t1126 + t1128) - 0.1904761905e0 * t1342 * t442 + t30 * (t1172 + t1174) - 0.1904761905e0 * t1347 * t442 + t12 * (t1231 + t1233) - 0.9523809525e-1 * t1352 * t442 + t1 * (t1292 + t1294) - 0.1904761905e-1 * t1357 * t442 + omega * t373 * t1099;
        double t1401 = 0.9523809525e-2 * t447 * t431 - 0.1904761905e-1 * t137 * omega * t439 + t148 * (0.4761904762e-1 * t462 + t1109 + t1110) - 0.9523809525e-1 * t3 * t1299 * t439 + t13 * (0.9523809525e-1 * t507 + t1146) - 0.1904761905e0 * t51 * t1238 * t439 + t30 * (0.9523809525e-1 * t571 + t1202) - 0.1904761905e0 * t76 * t1179 * t439 + t12 * (0.4761904762e-1 * t638 + t1265) - 0.9523809525e-1 * t100 * t1133 * t439 + t1 * (0.9523809525e-2 * t699 + t1328) - 0.1904761905e-1 * t373 * t1105 * t439 + t373 * (t1109 + t1110);
        double t1403 = 0.9523809525e-2 * t57 * t82 * t130 + 0.2857142858e-1 * t734 * t431 + t737 * (0.5714285715e-1 * t1 * t130 - 0.1904761905e-1 * x * t442 + t1097 + t1100) + t364 * (0.1523809524e0 * t1 * t431 - 0.1904761905e-1 * t1105 * t439 + 0.1428571429e0 * t462 + t1109 + t1110) + t377 * (0.1428571429e0 * t12 * t130 - 0.9523809525e-1 * t136 * t442 + t1 * t1121 - 0.9523809525e-1 * t1123 * t442 + t1126 + t1128) + t291 * (0.3333333334e0 * t12 * t431 - 0.9523809525e-1 * t1133 * t439 + t1 * (0.6190476191e0 * t462 + t1137 + t1138) - 0.9523809525e-1 * t3 * t1105 * t439 + 0.2857142858e0 * t507 + t1146) + t82 * (0.1904761905e0 * t30 * t130 - 0.1904761905e0 * t252 * t442 + t12 * t1157 - 0.3809523810e0 * t1159 * t442 + t1 * t1167 - 0.1904761905e0 * t1169 * t442 + t1172 + t1174) + t238 * (0.3809523810e0 * t30 * t431 - 0.1904761905e0 * t1179 * t439 + t12 * (0.1047619048e1 * t462 + t1183 + t1184) - 0.3809523810e0 * t3 * t1133 * t439 + t1 * (0.9523809525e0 * t507 + t1194) - 0.1904761905e0 * t51 * t1105 * t439 + 0.2857142858e0 * t571 + t1202) + t57 * (0.1428571429e0 * t13 * t130 - 0.1904761905e0 * t201 * t442 + t30 * t1157 - 0.5714285715e0 * t1210 * t442 + t12 * (omega * t3 * t1215 + t1213) - 0.5714285715e0 * t1220 * t442 + t1 * t1226 - 0.1904761905e0 * t1228 * t442 + t1231 + t1233) + t187 * (0.2380952381e0 * t13 * t431 - 0.1904761905e0 * t1238 * t439 + t30 * (0.8571428572e0 * t462 + t1183 + t1184) - 0.5714285715e0 * t3 * t1179 * t439 + t12 * (0.1142857143e1 * t507 + t3 * t1250) - 0.5714285715e0 * t51 * t1133 * t439 + t1 * (0.6666666668e0 * t571 + t1258) - 0.1904761905e0 * t76 * t1105 * t439 + 0.1428571429e0 * t638 + t1265) + t32 * t1295 + t150 * t1329 + t2 * t1362 + y * t1401;
        double t1411 = Vx * rho_0;
        double t1413 = rho_0 * omega;
        double t1414 = y * t1413;
        double t1416 = -0.6666666665e-1 * t1411 - 0.6666666665e-1 * t1414;
        double t1417 = t72 * t100;
        double t1421 = -0.3333333332e0 * t1414 - 0.3333333332e0 * t1411;
        double t1423 = t2 * t1411;
        double t1424 = 0.3333333332e0 * t1423;
        double t1425 = t150 * t1413;
        double t1426 = 0.3333333332e0 * t1425;
        double t1428 = z * t100;
        double t1432 = -0.6666666665e0 * t1414 - 0.6666666665e0 * t1411;
        double t1436 = -0.1333333333e1 * t1425 - 0.1333333333e1 * t1423;
        double t1438 = t32 * t1411;
        double t1439 = 0.6666666665e0 * t1438;
        double t1440 = t187 * t1413;
        double t1441 = 0.6666666665e0 * t1440;
        double t1443 = t72 * t51;
        double t1452 = -0.2000000000e1 * t1438 - 0.2000000000e1 * t1440;
        double t1454 = t57 * t1411;
        double t1455 = 0.6666666665e0 * t1454;
        double t1456 = t238 * t1413;
        double t1457 = 0.6666666665e0 * t1456;
        double t1468 = 0.3333333332e0 * t291 * t1413;
        double t1470 = 0.3333333332e0 * t82 * t1411;
        double t1495 = 0.1200000000e1 * t24;
        double t1496 = 0.3999999999e1 * t21;
        double t1506 = 0.9999999998e0 * t21;
        double t1509 = 0.1200000000e1 * t34;
        double t1511 = t2 * t41 * t40;
        double t1527 = 0.2999999999e1 * t21;
        double t1554 = t782 * t770;
        double t1557 = omega * Vx;
        double t1558 = y * t1557;
        double t1560 = 0.5e0 * t922;
        double t1561 = x * t143;
        double t1562 = 0.5e0 * t122;
        double t1563 = 0.5e0 * t123;
        double t1567 = t780 * t782 * t772;
        double t1569 = t772 * rho_0;
        double t1570 = t761 + Vx;
        double t1576 = 0.5000000000e0 * t782 * t780 * (-t1105 + Vy) * t1570 * t1569;
        double t1580 = t782 * Vz * t1570 * t780 * t1569;
        double t1586 = t150 * mu * t11;
        double t1592 = t187 * mu * t11;
        double t1598 = t238 * mu;
        double t1599 = t1598 * t11;
        double t1620 = t122 * H;
        double t1622 = 0.9523809525e-2 * t121;
        double t1623 = t1622 + t1046;
        double t1625 = 0.1904761905e-1 * t1620 + H * t1623;
        double t1631 = 0.4761904762e-1 * t121;
        double t1632 = t1631 + t868;
        double t1634 = 0.9523809525e-1 * t1620 + H * t1632;
        double t1636 = omega * t100 * t1634;
        double t1640 = 0.9523809525e-1 * t121;
        double t1641 = t852 + t1640;
        double t1643 = 0.1904761905e0 * t1620 + H * t1641;
        double t1645 = omega * t51 * t1643;
        double t1648 = 0.9523809525e-1 * t424 * t127 + 0.1904761905e0 * t412 * t127 + 0.1904761905e-1 * t390 * t127 + 0.1904761905e0 * t418 * t127 + 0.9523809525e-1 * t421 * t127 - 0.3809523810e-1 * t977 + 0.1904761905e-1 * t979 + omega * t373 * t1625 + 0.9523809525e-2 * t1078 + 0.1904761905e-1 * t1080 + t2 * (t972 + t1636) + t57 * (t848 + t1645);
        double t1650 = omega * t76 * t1643;
        double t1653 = omega * t1625;
        double t1657 = omega * t3 * t1634;
        double t1670 = t384 * t139;
        double t1678 = omega * t1634;
        double t1684 = 0.1904761905e0 * t121;
        double t1685 = t1684 + t881;
        double t1687 = 0.3809523810e0 * t1620 + H * t1685;
        double t1689 = omega * t3 * t1687;
        double t1695 = t1154 + t895;
        double t1697 = 0.5714285715e0 * t1620 + H * t1695;
        double t1699 = omega * t51 * t1697;
        double t1710 = 0.6666666668e-1 * t377 * t131 + 0.9523809525e-1 * t333 * t127 + t82 * (t865 + t1678) + 0.3809523810e0 * t321 * t127 + t57 * (t879 + t1689) + 0.5714285715e0 * t324 * t127 + t32 * (t893 + t1699) + 0.3809523810e0 * t327 * t127 + t2 * (omega * t76 * t1687 + t908) + 0.9523809525e-1 * t330 * t127 + t917 + t1636;
        double t1720 = t297 * t139;
        double t1728 = omega * t1643;
        double t1749 = t232 * t139;
        double t1771 = t196 * t139;
        double t1789 = t153 * t139;
        double t1799 = t32 * (t967 + t1650) + t377 * (t1072 + t1653) + t82 * (t1075 + t1657) + x * (-0.1904761905e0 * t358 * t139 - 0.3809523810e0 * t361 * t139 - 0.1904761905e0 * t355 * t139 - 0.3809523810e0 * t352 * t139 - 0.3809523810e-1 * t378 * t139 - 0.3809523810e-1 * t1670) + t1 * t1710 + t136 * (-0.1142857143e1 * t288 * t139 - 0.7619047620e0 * t294 * t139 - 0.7619047620e0 * t285 * t139 - 0.1904761905e0 * t303 * t139 - 0.1904761905e0 * t1720) + t12 * (0.1904761905e0 * t82 * t131 + 0.1904761905e0 * t269 * t127 + t57 * (t1005 + t1728) + 0.5714285715e0 * t266 * t127 + t32 * (omega * t3 * t1697 + t1011) + 0.5714285715e0 * t263 * t127 + t2 * (t1019 + t1699) + 0.1904761905e0 * t260 * t127 + t1025 + t1650) + t252 * (-0.1142857143e1 * t229 * t139 - 0.1142857143e1 * t226 * t139 - 0.3809523810e0 * t1749 - 0.3809523810e0 * t247 * t139) + t30 * (0.2857142858e0 * t57 * t131 + 0.1904761905e0 * t203 * t127 + t32 * (t985 + t1728) + 0.3809523810e0 * t212 * t127 + t2 * (t992 + t1689) + 0.1904761905e0 * t209 * t127 + t998 + t1645) + t201 * (-0.7619047620e0 * t190 * t139 - 0.3809523810e0 * t1771 - 0.3809523810e0 * t193 * t139) + t13 * (0.2380952381e0 * t32 * t131 + 0.9523809525e-1 * t175 * t127 + t2 * (t956 + t1678) + 0.9523809525e-1 * t169 * t127 + t962 + t1657) + t164 * (-0.1904761905e0 * t156 * t139 - 0.1904761905e0 * t1789) + t148 * (0.1047619048e0 * t2 * t131 + 0.1904761905e-1 * t140 * t127 + t1043 + t1653);
        double t1813 = -t1622 - 0.2857142858e-1 * t122 - t1046;
        double t1814 = omega * t1813;
        double t1821 = t3 * t442;
        double t1823 = Vy * t122;
        double t1824 = 0.9523809525e-2 * t1823;
        double t1825 = Vy * t1623;
        double t1833 = -t1631 - 0.1428571429e0 * t122 - t868;
        double t1835 = omega * t1833 - t1117;
        double t1840 = omega * t3 * t1833;
        double t1845 = t150 * omega;
        double t1849 = 0.4761904762e-1 * t1823;
        double t1850 = Vy * t1632;
        double t1856 = t51 * t442;
        double t1858 = t1849 + t1850;
        double t1859 = t3 * t1858;
        double t1866 = -t1640 - t852 - t896;
        double t1868 = omega * t1866 - t1153;
        double t1873 = -t1684 - 0.5714285715e0 * t122 - t881;
        double t1876 = omega * t3 * t1873 - t1162;
        double t1881 = omega * t51 * t1866;
        double t1886 = t187 * omega;
        double t1890 = 0.9523809525e-1 * t1823;
        double t1891 = Vy * t1641;
        double t1900 = 0.1904761905e0 * t1823 + Vy * t1685;
        double t1901 = t3 * t1900;
        double t1907 = t76 * t442;
        double t1909 = t1890 + t1891;
        double t1910 = t51 * t1909;
        double t1921 = -t1154 - 0.8571428572e0 * t122 - t895;
        double t1930 = omega * t51 * t1921 - t1223;
        double t1935 = omega * t76 * t1866;
        double t1940 = t238 * omega;
        double t1952 = 0.2857142858e0 * t1823 + Vy * t1695;
        double t1960 = t51 * t1952;
        double t1966 = t100 * t442;
        double t1968 = t76 * t1909;
        double t1991 = omega * t100 * t1833;
        double t1992 = -0.5714285715e-1 * t377 * t130 - 0.9523809525e-1 * t291 * t431 + t82 * t1835 - 0.3809523810e0 * t309 * t431 + t57 * t1876 - 0.5714285715e0 * t300 * t431 + t32 * t1930 - 0.3809523810e0 * t306 * t431 + t2 * (omega * t76 * t1873 - t1284) - 0.9523809525e-1 * t282 * t431 - t1292 + t1991;
        double t1996 = t291 * omega;
        double t2024 = t373 * t442;
        double t2026 = t100 * t1858;
        double t2027 = 0.7619047620e-1 * t377 * t442 + 0.9523809525e-1 * t1996 * t439 + t82 * (0.3333333334e0 * t1821 + t1849 + t1850) + 0.3809523810e0 * t3 * t1940 * t439 + t57 * (0.5714285715e0 * t1856 + t1901) + 0.5714285715e0 * t51 * t1886 * t439 + t32 * (0.4761904762e0 * t1907 + t1960) + 0.3809523810e0 * t76 * t1845 * t439 + t2 * (0.1904761905e0 * t1966 + t76 * t1900) + 0.9523809525e-1 * t100 * t761 * t439 + 0.2857142858e-1 * t2024 + t2026;
        double t2055 = -0.9523809525e-2 * t737 * t130 - 0.1904761905e-1 * t364 * t431 + t377 * (-t1097 + t1814) - 0.9523809525e-1 * t381 * t431 + t82 * (-t1126 + t1840) - 0.1904761905e0 * t367 * t431 + t57 * (-t1172 + t1881) - 0.1904761905e0 * t370 * t431 + t32 * (-t1231 + t1935) - 0.9523809525e-1 * t349 * t431 + t2 * (-t1292 + t1991) - 0.1904761905e-1 * t374 * t431 + omega * t373 * t1813;
        double t2094 = 0.9523809525e-2 * t737 * t442 + 0.1904761905e-1 * t364 * omega * t439 + t377 * (0.4761904762e-1 * t1821 + t1824 + t1825) + 0.9523809525e-1 * t3 * t1996 * t439 + t82 * (0.9523809525e-1 * t1856 + t1859) + 0.1904761905e0 * t51 * t1940 * t439 + t57 * (0.9523809525e-1 * t1907 + t1910) + 0.1904761905e0 * t76 * t1886 * t439 + t32 * (0.4761904762e-1 * t1966 + t1968) + 0.9523809525e-1 * t100 * t1845 * t439 + t2 * (0.9523809525e-2 * t2024 + t2026) + 0.1904761905e-1 * t373 * t761 * t439 + t373 * (t1824 + t1825);
        double t2096 = -0.9523809525e-2 * t30 * t13 * t130 + 0.2857142858e-1 * t437 * t442 + t447 * (-0.5714285715e-1 * t2 * t130 - 0.1904761905e-1 * y * t431 - t1097 + t1814) + t137 * (0.1523809524e0 * t2 * t442 + 0.1904761905e-1 * t761 * t439 + 0.1428571429e0 * t1821 + t1824 + t1825) + t148 * (-0.1428571429e0 * t32 * t130 - 0.9523809525e-1 * t150 * t431 + t2 * t1835 - 0.9523809525e-1 * t159 * t431 - t1126 + t1840) + t164 * (0.3333333334e0 * t32 * t442 + 0.9523809525e-1 * t1845 * t439 + t2 * (0.6190476191e0 * t1821 + t1849 + t1850) + 0.9523809525e-1 * t3 * t761 * t439 + 0.2857142858e0 * t1856 + t1859) + t13 * (-0.1904761905e0 * t57 * t130 - 0.1904761905e0 * t187 * t431 + t32 * t1868 - 0.3809523810e0 * t184 * t431 + t2 * t1876 - 0.1904761905e0 * t181 * t431 - t1172 + t1881) + t201 * (0.3809523810e0 * t57 * t442 + 0.1904761905e0 * t1886 * t439 + t32 * (0.1047619048e1 * t1821 + t1890 + t1891) + 0.3809523810e0 * t3 * t1845 * t439 + t2 * (0.9523809525e0 * t1856 + t1901) + 0.1904761905e0 * t51 * t761 * t439 + 0.2857142858e0 * t1907 + t1910) + t30 * (-0.1428571429e0 * t82 * t130 - 0.1904761905e0 * t238 * t431 + t57 * t1868 - 0.5714285715e0 * t241 * t431 + t32 * (omega * t3 * t1921 - t1213) - 0.5714285715e0 * t235 * t431 + t2 * t1930 - 0.1904761905e0 * t244 * t431 - t1231 + t1935) + t252 * (0.2380952381e0 * t82 * t442 + 0.1904761905e0 * t1940 * t439 + t57 * (0.8571428572e0 * t1821 + t1890 + t1891) + 0.5714285715e0 * t3 * t1886 * t439 + t32 * (0.1142857143e1 * t1856 + t3 * t1952) + 0.5714285715e0 * t51 * t1845 * t439 + t2 * (0.6666666668e0 * t1907 + t1960) + 0.1904761905e0 * t76 * t761 * t439 + 0.1428571429e0 * t1966 + t1968) + t12 * t1992 + t136 * t2027 + t1 * t2055 + x * t2094;
        double t2106 = t1 * mu;
        double t2107 = t2106 * t11;
        double t2112 = t12 * mu;
        double t2113 = t2112 * t11;
        double t2121 = 0.3999999999e0 * t53;
        double t2124 = t30 * mu;
        double t2125 = t2124 * t11;
        double t2158 = -0.2666666666e0 * t84 + t57 * (-0.7999999998e0 * t21 - 0.6666666665e0 * t2107 - 0.6666666665e0 * t27) + t32 * (-0.3999999999e0 * t2113 + t1 * (-0.7999999998e0 * t27 - 0.1400000000e1 * t21) + 0.4599999999e1 * t43 - 0.5999999998e1 * t49 - t2121) + t2 * (0.1333333333e0 * t2125 + t12 * (0.3999999999e0 * t27 - 0.3999999999e0 * t21) + t1 * (t2121 + 0.4199999999e1 * t43 - 0.4999999999e1 * t49) + 0.4599999999e1 * t70 + 0.1366666666e2 * t74 + 0.1333333333e0 * t78) + 0.1333333333e0 * t15 + t30 * (0.2000000000e0 * t21 + 0.5333333332e0 * t27) + t12 * (-0.3999999999e0 * t43 + 0.9999999998e0 * t49 + 0.7999999998e0 * t53) + t1 * (-0.1400000000e1 * t70 - 0.3333333332e0 * t74 + 0.5333333332e0 * t78) + 0.1333333333e0 * t102 - 0.7999999998e0 * t105 - 0.1333333333e1 * t109;
        double t2218 = t1 * omega;
        double t2219 = t3 * t2218;
        double t2255 = t51 * t2218;
        double t2258 = t12 * omega;
        double t2259 = t3 * t2258;
        double t2301 = t51 * t2258;
        double t2304 = t76 * t2218;
        double t2307 = t30 * omega;
        double t2308 = t3 * t2307;
        double t2337 = 0.6666666665e-1 * t1670 - 0.6666666665e0 * t1342 * t131 - 0.3333333332e0 * t1352 * t131 + 0.6666666665e0 * t936 * t139 + 0.3333333332e0 * t940 * t139 + 0.6666666665e0 * t932 * t139 + 0.3333333332e0 * t943 * t139 - 0.6666666665e-1 * t137 * t131 + 0.6666666665e-1 * t947 * t139 - 0.6666666665e-1 * t1357 * t131 - 0.3333333332e0 * t1337 * t131 - 0.6666666665e0 * t1347 * t131;
        double t2345 = t100 * t2218;
        double t2348 = t76 * t2258;
        double t2351 = t13 * omega;
        double t2352 = t3 * t2351;
        double t2355 = t51 * t2307;
        double t2369 = -0.6666666665e0 * t836 * t127 + 0.6666666665e-1 * t394 - 0.6666666665e-1 * t976 * t127 + 0.3333333332e0 * t2345 * t218 + 0.6666666665e0 * t2348 * t218 + 0.3333333332e0 * t2352 * t218 + 0.6666666665e0 * t2355 * t218 - 0.3333333332e0 * t832 * t127 - 0.3333333332e0 * t828 * t127 + 0.6666666665e-1 * t148 * t143 * t127 - 0.6666666665e0 * t840 * t127 - 0.6666666665e-1 * t844 * t127;
        double t2371 = t364 * rho_0 * (0.6666666665e-1 * t113 * t139 - 0.6666666665e-1 * x * t131) + t377 * rho_0 * (0.6666666665e-1 * t144 - 0.6666666665e-1 * t843 * t127) + t291 * rho_0 * (-0.3333333332e0 * t136 * t131 + 0.3333333332e0 * t922 * t139 + 0.3333333332e0 * t1789 - 0.3333333332e0 * t1123 * t131) + t82 * rho_0 * (0.3333333332e0 * t1 * t143 * t127 - 0.3333333332e0 * t959 * t127 + 0.3333333332e0 * t167 - 0.3333333332e0 * t827 * t127) + t238 * rho_0 * (0.1333333333e1 * t923 * t139 - 0.6666666665e0 * t252 * t131 - 0.6666666665e0 * t1169 * t131 + 0.6666666665e0 * t1771 + 0.6666666665e0 * t927 * t139 - 0.1333333333e1 * t1159 * t131) + t57 * rho_0 * (-0.6666666665e0 * t995 * t127 + 0.6666666665e0 * t12 * t143 * t127 - 0.1333333333e1 * t989 * t127 + 0.6666666665e0 * t216 - 0.6666666665e0 * t835 * t127 + 0.1333333333e1 * t2219 * t218) + t187 * rho_0 * (0.6666666665e0 * t935 * t139 + 0.2000000000e1 * t1031 * t139 + 0.2000000000e1 * t1028 * t139 - 0.6666666665e0 * t201 * t131 - 0.6666666665e0 * t1228 * t131 - 0.2000000000e1 * t1220 * t131 - 0.2000000000e1 * t1210 * t131 + 0.6666666665e0 * t1749) + t32 * rho_0 * (-0.2000000000e1 * t1008 * t127 + 0.6666666665e0 * t30 * t143 * t127 + 0.6666666665e0 * t255 - 0.6666666665e0 * t1022 * t127 - 0.2000000000e1 * t1016 * t127 - 0.6666666665e0 * t839 * t127 + 0.2000000000e1 * t2255 * t218 + 0.2000000000e1 * t2259 * t218) + t150 * rho_0 * (-0.3333333332e0 * t1289 * t131 - 0.1333333333e1 * t1273 * t131 - 0.2000000000e1 * t1277 * t131 + 0.1333333333e1 * t1061 * t139 + 0.1333333333e1 * t1058 * t139 + 0.2000000000e1 * t1064 * t139 - 0.3333333332e0 * t164 * t131 + 0.3333333332e0 * t939 * t139 + 0.3333333332e0 * t1720 - 0.1333333333e1 * t1281 * t131) + t2 * rho_0 * (-0.2000000000e1 * t890 * t127 - 0.3333333332e0 * t913 * t127 + 0.3333333332e0 * t319 + 0.3333333332e0 * t13 * t143 * t127 - 0.1333333333e1 * t904 * t127 - 0.1333333333e1 * t875 * t127 - 0.3333333332e0 * t831 * t127 + 0.2000000000e1 * t2301 * t218 + 0.1333333333e1 * t2304 * t218 + 0.1333333333e1 * t2308 * t218) + y * rho_0 * t2337 + rho_0 * t2369;
        double t2390 = 0.6666666665e-1 * t121;
        double t2392 = 0.6666666665e-1 * t123;
        double t2393 = t2390 + 0.2000000000e0 * t122 + t2392;
        double t2395 = 0.3333333332e0 * t453 + omega * t2393;
        double t2398 = 0.6666666665e-1 * t1823;
        double t2400 = Vy * (-t2392 - t2390);
        double t2419 = 0.3333333332e0 * t123;
        double t2420 = 0.3333333332e0 * t121;
        double t2422 = t2419 + t2420 + 0.9999999998e0 * t122;
        double t2424 = 0.1666666666e1 * t453 + omega * t2422;
        double t2427 = 0.3333333332e0 * t1823;
        double t2429 = Vy * (-t2420 - t2419);
        double t2435 = 0.6666666665e0 * t501 + omega * t3 * t2422;
        double t2438 = -t2427 + t2429;
        double t2439 = t3 * t2438;
        double t2462 = 0.2000000000e1 * t122;
        double t2463 = 0.6666666665e0 * t123;
        double t2464 = 0.6666666665e0 * t121;
        double t2465 = t2462 + t2463 + t2464;
        double t2467 = 0.3333333332e1 * t453 + omega * t2465;
        double t2470 = 0.6666666665e0 * t1823;
        double t2472 = Vy * (-t2464 - t2463);
        double t2476 = 0.1333333333e1 * t123;
        double t2478 = 0.1333333333e1 * t121;
        double t2479 = t2476 + 0.3999999999e1 * t122 + t2478;
        double t2482 = 0.2666666666e1 * t501 + omega * t3 * t2479;
        double t2488 = -0.1333333333e1 * t1823 + Vy * (-t2478 - t2476);
        double t2489 = t3 * t2488;
        double t2495 = 0.6666666665e0 * t565 + omega * t51 * t2465;
        double t2498 = -t2470 + t2472;
        double t2499 = t51 * t2498;
        double t2530 = 0.2000000000e1 * t121;
        double t2532 = 0.2000000000e1 * t123;
        double t2533 = t2530 + 0.5999999998e1 * t122 + t2532;
        double t2542 = -0.2000000000e1 * t1823 + Vy * (-t2530 - t2532);
        double t2549 = 0.2000000000e1 * t565 + omega * t51 * t2533;
        double t2552 = t51 * t2542;
        double t2558 = 0.3333333332e0 * t632 + omega * t76 * t2465;
        double t2561 = t76 * t2498;
        double t2562 = 0.9999999998e0 * t164 * t130 - 0.2333333333e1 * t13 * t442 + t201 * t2467 + t30 * (-0.7333333332e1 * t1821 - t2470 + t2472) + t252 * (0.3999999999e1 * t501 + omega * t3 * t2533) + t12 * (-0.7999999998e1 * t1856 + t3 * t2542) + t136 * t2549 + t1 * (-0.3333333332e1 * t1907 + t2552) + x * t2558 - 0.3333333332e0 * t1966 + t2561;
        double t2615 = 0.6666666665e-1 * t693 + omega * t100 * t2422;
        double t2618 = t100 * t2438;
        double t2619 = 0.3999999999e0 * t137 * t130 - 0.1066666666e1 * t148 * t442 + t164 * t2424 + t13 * (-0.4333333332e1 * t1821 - t2427 + t2429) + t201 * t2482 + t30 * (-0.6666666665e1 * t1856 + t2489) + t252 * t2549 + t12 * (-0.4666666666e1 * t1907 + t2552) + t136 * (0.6666666665e0 * t632 + omega * t76 * t2479) + t1 * (-0.1333333333e1 * t1966 + t76 * t2488) + x * t2615 - 0.6666666665e-1 * t2024 + t2618;
        double t2628 = t148 * omega;
        double t2646 = -0.6666666665e0 * t2352 * t439 - 0.1333333333e1 * t2348 * t439 - 0.6666666665e0 * t2345 * t439 - 0.1333333333e0 * t2628 * t439 + 0.1333333333e1 * t1347 * t431 + 0.6666666665e0 * t1337 * t431 + 0.1333333333e0 * t1357 * t431 + 0.6666666665e0 * t1352 * t431 - 0.1333333333e0 * t718 + 0.1333333333e1 * t1342 * t431 - 0.1333333333e1 * t2355 * t439 + 0.1333333333e0 * t137 * t431;
        double t2677 = 0.6666666665e-1 * t437 * t130 - 0.2000000000e0 * t447 * t442 + t137 * t2395 + t148 * (-0.9999999998e0 * t1821 - t2398 + t2400) + t164 * t2435 + t13 * (-0.2000000000e1 * t1856 + t2439) + t201 * t2495 + t30 * (-0.2000000000e1 * t1907 + t2499) + t252 * t2558 + t12 * (-0.9999999998e0 * t1966 + t2561) + t136 * t2615 + t1 * (-0.2000000000e0 * t2024 + t2618) + t1105 * t373 * t2393 + t373 * (-t2398 + t2400);
        double t2680 = t8 * t2371 + t734 * rho_0 * (-0.6666666665e-1 * t442 + 0.6666666665e-1 * x * t130) + t737 * rho_0 * (0.1333333333e0 * x * t431 - 0.1333333333e0 * t440) + t364 * rho_0 * (0.3999999999e0 * t136 * t130 - 0.5333333332e0 * t1 * t442 + x * t2395 - 0.3333333332e0 * t1821 - t2398 + t2400) + t377 * rho_0 * (0.6666666665e0 * t1123 * t431 - 0.6666666665e0 * t474 - 0.6666666665e0 * t2218 * t439 + 0.6666666665e0 * t136 * t431) + t291 * rho_0 * (0.9999999998e0 * t252 * t130 - 0.1666666666e1 * t12 * t442 + t136 * t2424 + t1 * (-0.2333333333e1 * t1821 - t2427 + t2429) + x * t2435 - 0.6666666665e0 * t1856 + t2439) + t82 * rho_0 * (-0.2666666666e1 * t2219 * t439 + 0.2666666666e1 * t1159 * t431 - 0.1333333333e1 * t525 - 0.1333333333e1 * t2258 * t439 + 0.1333333333e1 * t1169 * t431 + 0.1333333333e1 * t252 * t431) + t238 * rho_0 * (0.1333333333e1 * t201 * t130 - 0.2666666666e1 * t30 * t442 + t252 * t2467 + t12 * (-0.5999999998e1 * t1821 - t2470 + t2472) + t136 * t2482 + t1 * (-0.3999999999e1 * t1856 + t2489) + x * t2495 - 0.6666666665e0 * t1907 + t2499) + t57 * rho_0 * (-0.3999999999e1 * t2255 * t439 + 0.3999999999e1 * t1210 * t431 + 0.1333333333e1 * t201 * t431 - 0.1333333333e1 * t2307 * t439 + 0.1333333333e1 * t1228 * t431 + 0.3999999999e1 * t1220 * t431 - 0.1333333333e1 * t587 - 0.3999999999e1 * t2259 * t439) + t187 * rho_0 * t2562 + t32 * rho_0 * (-0.2666666666e1 * t2308 * t439 - 0.3999999999e1 * t2301 * t439 - 0.2666666666e1 * t2304 * t439 + 0.2666666666e1 * t1281 * t431 + 0.2666666666e1 * t1273 * t431 + 0.6666666665e0 * t164 * t431 + 0.6666666665e0 * t1289 * t431 - 0.6666666665e0 * t2351 * t439 + 0.3999999999e1 * t1277 * t431 - 0.6666666665e0 * t657) + t150 * rho_0 * t2619 + t2 * rho_0 * t2646 + y * rho_0 * t2677;
        double t2686 = Vy * rho_0;
        double t2688 = x * t1413;
        double t2690 = -0.6666666665e-1 * t2686 + 0.6666666665e-1 * t2688;
        double t2694 = -0.3333333332e0 * t2686 + 0.3333333332e0 * t2688;
        double t2696 = t136 * t1413;
        double t2697 = 0.3333333332e0 * t2696;
        double t2698 = t1 * t2686;
        double t2699 = 0.3333333332e0 * t2698;
        double t2704 = -0.6666666665e0 * t2686 + 0.6666666665e0 * t2688;
        double t2708 = -0.1333333333e1 * t2698 + 0.1333333333e1 * t2696;
        double t2710 = t252 * t1413;
        double t2711 = 0.6666666665e0 * t2710;
        double t2712 = t12 * t2686;
        double t2713 = 0.6666666665e0 * t2712;
        double t2723 = -0.2000000000e1 * t2712 + 0.2000000000e1 * t2710;
        double t2725 = t201 * t1413;
        double t2726 = 0.6666666665e0 * t2725;
        double t2727 = t30 * t2686;
        double t2728 = 0.6666666665e0 * t2727;
        double t2739 = 0.3333333332e0 * t164 * t1413;
        double t2741 = 0.3333333332e0 * t13 * t2686;
        double t2767 = 0.1200000000e1 * t2107;
        double t2780 = t1 * t41 * t40;
        double t2782 = 0.1200000000e1 * t2113;
        double t2826 = 0.5e0 * t156;
        double t2827 = 0.5e0 * t121;
        double t2838 = t782 * t780 * Vz * t772 * rho_0 * (0.5e0 * t1105 - 0.5e0 * Vy);
        double t2839 = t780 * t10;
        double t2846 = ball * t46;
        double t2848 = mu * t2846 * t45;
        double t2849 = 0.9999999998e0 * t2848;
        double t2851 = t72 * mu * t11;
        double t2852 = 0.1200000000e1 * t2851;
        double t2853 = t48 * t40;
        double t2854 = 0.2999999999e1 * t2853;
        double t2863 = -0.5999999998e1 * t2853 + 0.2000000000e1 * t2848 - 0.2399999999e1 * t2851;
        double t2867 = 0.1200000000e1 * t107 * mu * t11;
        double t2869 = 0.1200000000e2 * t42 * t47;
        double t2871 = 0.9999999998e0 * t73 * t40;
        double t2877 = -t2852 - t2854 + t2849;
        double t2879 = -t2867 - t2869 + t2871;
        double t2883 = 0.3999999999e0 * t1443 * mu * t11;
        double t2885 = 0.3999999999e1 * t108 * t40;
        double t2887 = 0.7999999998e1 * t69 * t47;
        double t2892 = H * Vz;
        double t2893 = rho_0 * t113;
        double t2901 = t113 * H;
        double t2903 = 0.3333333332e0 * t3 * t2901;
        double t2916 = t13 * rho_0;
        double t2921 = 0.6666666665e0 * t51 * t2901;
        double t2922 = t3 * t2;
        double t2937 = t30 * rho_0;
        double t2939 = t3 * t32;
        double t2945 = 0.6666666665e0 * t76 * t2901;
        double t2946 = t51 * t2;
        double t2963 = t12 * rho_0;
        double t2965 = t3 * t57;
        double t2971 = 0.3333333332e0 * t100 * t2901;
        double t2972 = t76 * t2;
        double t2975 = t51 * t32;
        double t2994 = t1 * rho_0;
        double t2999 = 0.6666666665e-1 * t373 * t2901;
        double t3000 = t76 * t32;
        double t3003 = t100 * t2;
        double t3006 = t3 * t82;
        double t3009 = t51 * t57;
        double t3031 = 0.6666666665e-1 * t137 * t2893 * t2892 - 0.6666666665e-1 * t148 * t1413 * Vz * t139 + t164 * rho_0 * Vz * (t2903 + 0.3333333332e0 * t2 * t2901) + t2916 * Vz * (-0.3333333332e0 * t219 * t139 - 0.3333333332e0 * t473 * t139) + t201 * rho_0 * Vz * (0.6666666665e0 * t32 * t2901 + t2921 + 0.1333333333e1 * t2922 * t2901) + t2937 * Vz * (-0.1333333333e1 * t220 * t139 - 0.6666666665e0 * t275 * t139 - 0.6666666665e0 * t524 * t139) + t252 * rho_0 * Vz * (0.2000000000e1 * t2939 * t2901 + 0.6666666665e0 * t57 * t2901 + t2945 + 0.2000000000e1 * t2946 * t2901) + t2963 * Vz * (-0.2000000000e1 * t276 * t139 - 0.2000000000e1 * t272 * t139 - 0.6666666665e0 * t342 * t139 - 0.6666666665e0 * t586 * t139) + t136 * rho_0 * Vz * (0.1333333333e1 * t2965 * t2901 + 0.3333333332e0 * t82 * t2901 + t2971 + 0.1333333333e1 * t2972 * t2901 + 0.2000000000e1 * t2975 * t2901) + t2994 * Vz * (-0.1333333333e1 * t339 * t139 - 0.1333333333e1 * t343 * t139 - 0.2000000000e1 * t336 * t139 - 0.3333333332e0 * t405 * t139 - 0.3333333332e0 * t656 * t139) + rho_0 * x * Vz * (0.6666666665e-1 * t377 * t2901 + t2999 + 0.6666666665e0 * t3000 * t2901 + 0.3333333332e0 * t3003 * t2901 + 0.3333333332e0 * t3006 * t2901 + 0.6666666665e0 * t3009 * t2901) + rho_0 * Vz * (-0.6666666665e-1 * t722 * t139 - 0.6666666665e-1 * t717 * t139 - 0.3333333332e0 * t406 * t139 - 0.6666666665e0 * t399 * t139 - 0.6666666665e0 * t396 * t139 - 0.3333333332e0 * t402 * t139);
        double t3033 = t113 * Vz;
        double t3041 = Vz * t123;
        double t3042 = 0.6666666665e-1 * t3041;
        double t3043 = 0.6666666665e-1 * t122;
        double t3046 = 0.3333333332e0 * t153;
        double t3059 = 0.3333333332e0 * t3;
        double t3066 = 0.3333333332e0 * t122;
        double t3067 = -0.1666666666e1 * t153 - t2420 - t3066;
        double t3071 = 0.6666666665e0 * t196;
        double t3072 = -t2420 - t3066;
        double t3073 = t3 * t3072;
        double t3087 = 0.6666666665e0 * t51;
        double t3095 = 0.6666666665e0 * t122;
        double t3096 = -0.3333333332e1 * t153 - t2464 - t3095;
        double t3102 = -t2478 - 0.1333333333e1 * t122;
        double t3104 = -0.2666666666e1 * t196 + t3 * t3102;
        double t3108 = 0.6666666665e0 * t232;
        double t3109 = -t2464 - t3095;
        double t3110 = t51 * t3109;
        double t3126 = 0.6666666665e0 * t76;
        double t3138 = -t2462 - t2530;
        double t3140 = -0.3999999999e1 * t196 + t3 * t3138;
        double t3146 = -0.2000000000e1 * t232 + t51 * t3138;
        double t3150 = 0.3333333332e0 * t297;
        double t3151 = t76 * t3109;
        double t3168 = 0.3333333332e0 * t100;
        double t3189 = -0.6666666665e0 * t297 + t76 * t3102;
        double t3193 = 0.6666666665e-1 * t384;
        double t3194 = t100 * t3072;
        double t3195 = -0.3999999999e0 * t378 - 0.6666666665e0 * t291 * t1557 + t82 * t3067 - 0.2666666666e1 * t309 * t1557 + t57 * t3104 - 0.3999999999e1 * t300 * t1557 + t32 * t3146 - 0.2666666666e1 * t306 * t1557 + t2 * t3189 - 0.6666666665e0 * t282 * t1557 - t3193 + t3194;
        double t3213 = 0.6666666665e-1 * t373;
        double t3225 = -t3046 - t2390 - t3043;
        double t3229 = -t3071 + t3073;
        double t3233 = -t3108 + t3110;
        double t3237 = -t3150 + t3151;
        double t3241 = -t3193 + t3194;
        double t3246 = t373 * (-t2390 - t3043);
        double t3247 = -0.6666666665e-1 * t737 * t113 - 0.1333333333e0 * t364 * t1557 + t377 * t3225 - 0.6666666665e0 * t381 * t1557 + t82 * t3229 - 0.1333333333e1 * t367 * t1557 + t57 * t3233 - 0.1333333333e1 * t370 * t1557 + t32 * t3237 - 0.6666666665e0 * t349 * t1557 + t2 * t3241 - 0.1333333333e0 * t374 * t1557 + t3246;
        double t3252 = t8 * t3031 - 0.6666666665e-1 * t437 * rho_0 * t3033 + 0.1333333333e0 * t447 * t1413 * Vy * Vz + t137 * rho_0 * (-t3042 + Vz * (-t3043 - t2390 - 0.1333333333e0 * t1558 - 0.3999999999e0 * t156 - t3046)) + t148 * rho_0 * Vz * (0.6666666665e0 * t172 + 0.6666666665e0 * t166) + t164 * rho_0 * (t3041 * (-0.3333333332e0 * t2 - t3059) + Vz * (-0.9999999998e0 * t193 - 0.6666666665e0 * t150 * t1557 + t2 * t3067 - 0.6666666665e0 * t159 * t1557 - t3071 + t3073)) + t2916 * Vz * (0.2666666666e1 * t2922 * t143 + 0.1333333333e1 * t206 + 0.1333333333e1 * t215) + t201 * rho_0 * (t3041 * (-0.6666666665e0 * t32 - t3087 - 0.1333333333e1 * t2922) + Vz * (-0.1333333333e1 * t247 - 0.1333333333e1 * t187 * t1557 + t32 * t3096 - 0.2666666666e1 * t184 * t1557 + t2 * t3104 - 0.1333333333e1 * t181 * t1557 - t3108 + t3110)) + t2937 * Vz * (0.3999999999e1 * t2939 * t143 + 0.3999999999e1 * t2946 * t143 + 0.1333333333e1 * t257 + 0.1333333333e1 * t254) + t252 * rho_0 * (t3041 * (-0.6666666665e0 * t57 - t3126 - 0.2000000000e1 * t2939 - 0.2000000000e1 * t2946) + Vz * (-0.9999999998e0 * t303 - 0.1333333333e1 * t238 * t1557 + t57 * t3096 - 0.3999999999e1 * t241 * t1557 + t32 * t3140 - 0.3999999999e1 * t235 * t1557 + t2 * t3146 - 0.1333333333e1 * t244 * t1557 - t3150 + t3151)) + t2963 * Vz * (0.2666666666e1 * t2965 * t143 + 0.3999999999e1 * t2975 * t143 + 0.2666666666e1 * t2972 * t143 + 0.6666666665e0 * t318 + 0.6666666665e0 * t315) + t136 * rho_0 * (t3041 * (-t3168 - 0.3333333332e0 * t82 - 0.2000000000e1 * t2975 - 0.1333333333e1 * t2972 - 0.1333333333e1 * t2965) + Vz * t3195) + t2994 * Vz * (0.6666666665e0 * t3003 * t143 + 0.1333333333e0 * t415 + 0.1333333333e1 * t3009 * t143 + 0.1333333333e1 * t3000 * t143 + 0.6666666665e0 * t3006 * t143 + 0.1333333333e0 * t393) + x * rho_0 * (t3041 * (-t3213 - 0.6666666665e-1 * t377 - 0.3333333332e0 * t3006 - 0.6666666665e0 * t3009 - 0.6666666665e0 * t3000 - 0.3333333332e0 * t3003) + Vz * t3247);
        double t3301 = t82 * rho_0;
        double t3305 = t3 * t1;
        double t3320 = t57 * rho_0;
        double t3324 = t51 * t1;
        double t3327 = t3 * t12;
        double t3344 = t32 * rho_0;
        double t3348 = t3 * t30;
        double t3351 = t76 * t1;
        double t3354 = t51 * t12;
        double t3373 = t2 * rho_0;
        double t3377 = t100 * t1;
        double t3380 = t51 * t30;
        double t3383 = t76 * t12;
        double t3386 = t3 * t13;
        double t3408 = 0.6666666665e-1 * t364 * t2893 * t2892 + 0.6666666665e-1 * t377 * t1413 * Vz * t127 + t291 * rho_0 * Vz * (t2903 + 0.3333333332e0 * t1 * t2901) + t3301 * Vz * (0.3333333332e0 * t2218 * t127 + 0.3333333332e0 * t473 * t127) + t238 * rho_0 * Vz * (0.6666666665e0 * t12 * t2901 + t2921 + 0.1333333333e1 * t3305 * t2901) + t3320 * Vz * (0.1333333333e1 * t2219 * t127 + 0.6666666665e0 * t2258 * t127 + 0.6666666665e0 * t524 * t127) + t187 * rho_0 * Vz * (t2945 + 0.6666666665e0 * t30 * t2901 + 0.2000000000e1 * t3324 * t2901 + 0.2000000000e1 * t3327 * t2901) + t3344 * Vz * (0.2000000000e1 * t2255 * t127 + 0.2000000000e1 * t2259 * t127 + 0.6666666665e0 * t2307 * t127 + 0.6666666665e0 * t586 * t127) + t150 * rho_0 * Vz * (0.3333333332e0 * t13 * t2901 + t2971 + 0.1333333333e1 * t3348 * t2901 + 0.1333333333e1 * t3351 * t2901 + 0.2000000000e1 * t3354 * t2901) + t3373 * Vz * (0.1333333333e1 * t2304 * t127 + 0.1333333333e1 * t2308 * t127 + 0.2000000000e1 * t2301 * t127 + 0.3333333332e0 * t2351 * t127 + 0.3333333332e0 * t656 * t127) + rho_0 * y * Vz * (t2999 + 0.6666666665e-1 * t148 * t2901 + 0.3333333332e0 * t3377 * t2901 + 0.6666666665e0 * t3380 * t2901 + 0.6666666665e0 * t3383 * t2901 + 0.3333333332e0 * t3386 * t2901) + rho_0 * Vz * (0.6666666665e-1 * t2628 * t127 + 0.6666666665e-1 * t717 * t127 + 0.3333333332e0 * t2352 * t127 + 0.3333333332e0 * t2345 * t127 + 0.6666666665e0 * t2355 * t127 + 0.6666666665e0 * t2348 * t127);
        double t3540 = -0.3999999999e0 * t947 + 0.6666666665e0 * t164 * t143 + t13 * t3067 + 0.2666666666e1 * t1273 * t143 + t30 * t3104 + 0.3999999999e1 * t1277 * t143 + t12 * t3146 + 0.2666666666e1 * t1281 * t143 + t1 * t3189 + 0.6666666665e0 * t1289 * t143 - t3193 + t3194;
        double t3586 = -0.6666666665e-1 * t447 * t113 + 0.1333333333e0 * t137 * t143 + t148 * t3225 + 0.6666666665e0 * t1337 * t143 + t13 * t3229 + 0.1333333333e1 * t1342 * t143 + t30 * t3233 + 0.1333333333e1 * t1347 * t143 + t12 * t3237 + 0.6666666665e0 * t1352 * t143 + t1 * t3241 + 0.1333333333e0 * t1357 * t143 + t3246;
        double t3591 = t8 * t3408 - 0.6666666665e-1 * t734 * rho_0 * t3033 - 0.1333333333e0 * t737 * t1413 * Vx * Vz + t364 * rho_0 * (-t3042 + Vz * (-0.3999999999e0 * t922 + 0.1333333333e0 * t1561 - t3046 - t2390 - t3043)) + t377 * rho_0 * Vz * (-0.6666666665e0 * t1 * t1557 - 0.6666666665e0 * t3 * t1557) + t291 * rho_0 * (t3041 * (-t3059 - 0.3333333332e0 * t1) + Vz * (-0.9999999998e0 * t927 + 0.6666666665e0 * t136 * t143 + t1 * t3067 + 0.6666666665e0 * t1123 * t143 - t3071 + t3073)) + t3301 * Vz * (-0.2666666666e1 * t3305 * t1557 - 0.1333333333e1 * t12 * t1557 - 0.1333333333e1 * t51 * t1557) + t238 * rho_0 * (t3041 * (-0.1333333333e1 * t3305 - t3087 - 0.6666666665e0 * t12) + Vz * (-0.1333333333e1 * t935 + 0.1333333333e1 * t252 * t143 + t12 * t3096 + 0.2666666666e1 * t1159 * t143 + t1 * t3104 + 0.1333333333e1 * t1169 * t143 - t3108 + t3110)) + t3320 * Vz * (-0.1333333333e1 * t30 * t1557 - 0.1333333333e1 * t76 * t1557 - 0.3999999999e1 * t3324 * t1557 - 0.3999999999e1 * t3327 * t1557) + t187 * rho_0 * (t3041 * (-0.2000000000e1 * t3327 - 0.2000000000e1 * t3324 - t3126 - 0.6666666665e0 * t30) + Vz * (-0.9999999998e0 * t939 + 0.1333333333e1 * t201 * t143 + t30 * t3096 + 0.3999999999e1 * t1210 * t143 + t12 * t3140 + 0.3999999999e1 * t1220 * t143 + t1 * t3146 + 0.1333333333e1 * t1228 * t143 - t3150 + t3151)) + t3344 * Vz * (-0.6666666665e0 * t13 * t1557 - 0.6666666665e0 * t100 * t1557 - 0.2666666666e1 * t3348 * t1557 - 0.2666666666e1 * t3351 * t1557 - 0.3999999999e1 * t3354 * t1557) + t150 * rho_0 * (t3041 * (-0.1333333333e1 * t3348 - 0.2000000000e1 * t3354 - 0.1333333333e1 * t3351 - 0.3333333332e0 * t13 - t3168) + Vz * t3540) + t3373 * Vz * (-0.6666666665e0 * t3386 * t1557 - 0.1333333333e0 * t373 * t1557 - 0.1333333333e0 * t148 * t1557 - 0.1333333333e1 * t3383 * t1557 - 0.1333333333e1 * t3380 * t1557 - 0.6666666665e0 * t3377 * t1557) + y * rho_0 * (t3041 * (-t3213 - 0.6666666665e-1 * t148 - 0.3333333332e0 * t3386 - 0.6666666665e0 * t3380 - 0.6666666665e0 * t3383 - 0.3333333332e0 * t3377) + Vz * t3586);
        double t3598 = Vz * t772;
        double t3599 = rho_0 * t3598;
        double t3612 = t125 * t30;
        double t3621 = z * t2;
        double t3625 = t125 * t12;
        double t3634 = mu * t11;
        double t3635 = t2 * t12;
        double t3639 = t32 * t1;
        double t3653 = 0.5e1 * t107 * t57 * t125 * t3599 + 0.5e1 * t1443 * t32 * t125 * t3599 + 0.25e1 * t1428 * t2 * t125 * t3599 + 0.5e1 * t107 * t3612 * t3599 + 0.63e2 * t2922 * t2106 * t19 * t17 - 0.75e2 * t3621 * t2106 * t2846 * t45 + 0.5e1 * z * t57 * t3625 * t3599 + 0.15e2 * t72 * t32 * t3625 * t3599 - 0.3e1 * t3 * t3635 * t3634 - 0.3e1 * t3 * t3639 * t3634 + 0.6e1 * t51 * t2 * t1 * t3634 - 0.45e1 * t30 * t41 * t40 - 0.45e1 * t57 * t41 * t40;
        double t3661 = t125 * t13;
        double t3665 = t125 * t1;
        double t3686 = 0.5e1 * t1443 * t3625 * t3599 + 0.5e0 * z * t125 * t148 * t3599 + 0.25e1 * t72 * t3661 * t3599 + 0.25e1 * t1428 * t3665 * t3599 + 0.5e0 * z * t377 * t125 * t3599 + 0.25e1 * t72 * t82 * t125 * t3599 - t84 + 0.2e1 * t102 - t15 - 0.12e2 * t105 - 0.20e2 * t109 - 0.135e2 * t3635 * t41 * t40 + 0.315e2 * t3327 * t41 * t40;
        double t3708 = z * t32;
        double t3712 = t72 * t2;
        double t3735 = -0.135e2 * t3639 * t41 * t40 + 0.24e2 * t3324 * t41 * t40 + 0.24e2 * t2946 * t41 * t40 + 0.315e2 * t2939 * t41 * t40 - 0.375e2 * z * t12 * t41 * t47 + 0.100e3 * t72 * t1 * t41 * t47 - 0.375e2 * t3708 * t41 * t47 + 0.100e3 * t3712 * t41 * t47 + 0.15e2 * t107 * t32 * t3665 * t3599 + 0.10e2 * t1443 * t2 * t3665 * t3599 + 0.5e1 * t3708 * t3612 * t3599 + 0.15e2 * t107 * t2 * t3625 * t3599 + 0.25e1 * z * t82 * t3665 * t3599;
        double t3775 = 0.10e2 * t72 * t57 * t3665 * t3599 + 0.25e1 * t3621 * t3661 * t3599 + 0.10e2 * t3712 * t3612 * t3599 - 0.4e1 * t2 * t2124 * t11 - t3 * t2124 * t11 - 0.6e1 * t32 * t2112 * t11 + 0.3e1 * t51 * t2112 * t11 - 0.4e1 * t57 * t2106 * t11 + 0.5e1 * t76 * t2106 * t11 - t3 * t58 * t11 + 0.3e1 * t51 * t33 * t11 + 0.5e1 * t76 * t23 * t11 + 0.5e0 * t1417 * t125 * rho_0 * t3598;
        A(0,3) = 1;
        A(1,4) = 1;
        A(2,5) = 1;
        A(3,0) = (0.4500000000e2 * t783 * (t125 * t111 + t772 * t768) * t10);
        A(3,1) = (0.5250000000e2 * t783 * (t125 * (0.5714285715e-1 * t238 * t786 * t11 + t187 * (0.1714285714e0 * t791 + x * t795) + t150 * (0.1714285714e0 * t800 + t136 * t804 + x * t809) + y * (0.5714285715e-1 * t814 + t252 * t795 + t136 * t809 + x * t821)) + t772 * (t8 * rho_0 * (t975 + t1084) + rho_0 * t1403)) * t10);
        A(3,2) = - (0.7499999998e1 * t1554 * t10 * (t779 * t772 * (t1417 * t1416 + t1428 * (t1 * t1421 - t1424 - t1426) + t1443 * (t1 * t1436 + t12 * t1432 - t1439 - t1441) + t107 * (t30 * t1432 + t12 * (-0.2000000000e1 * t1423 - 0.2000000000e1 * t1425) + t1 * t1452 - t1455 - t1457) + t72 * (t13 * t1421 + t30 * t1436 + t12 * t1452 + t1 * (-0.1333333333e1 * t1454 - 0.1333333333e1 * t1456) - t1468 - t1470) + z * (t148 * t1416 + t13 * (-t1424 - t1426) + t30 * (-t1439 - t1441) + t12 * (-t1455 - t1457) + t1 * (-t1468 - t1470) - 0.6666666665e-1 * t377 * t1411 - 0.6666666665e-1 * t364 * t1413)) - 0.3999999999e0 * t1443 * t786 * t11 + t107 * (-0.1200000000e1 * t791 + x * (-t1495 + t1496)) + 0.7999999998e1 * t1169 * t41 * t47 + t72 * (-0.1200000000e1 * t800 + t136 * (-0.2399999999e1 * t24 + t1506) + x * (-t1509 + 0.9999999998e0 * t1511)) + t3 * (-0.1200000000e2 * t2 * x * t41 * t47 - 0.1200000000e2 * t136 * t41 * t47) + z * (-0.3999999999e0 * t814 + t252 * (-t1495 - t1527) + t136 * (-t1509 - 0.5999999998e1 * t1511) + x * (-0.3999999999e0 * t59 - 0.2999999999e1 * t32 * t41 * t40)) + 0.9999999998e0 * t252 * t41 * t47 + 0.2000000000e1 * t2 * t136 * t41 * t47 + 0.9999999998e0 * t32 * x * t41 * t47));
        A(3,3) = - (t1567 * rho_0 * (t156 + 0.2e1 * t1558 + t121 + t1560 - t1561 + t1562 + t1563));
        A(3,4) = - t1576;
        A(3,5) = - (0.5e0 * t1580);
        A(4,0) = (0.5250000000e2 * t783 * t10 * (t125 * (0.5714285715e-1 * y * t813 * t11 + t252 * (0.1714285714e0 * t1586 + y * t795) + t136 * (0.1714285714e0 * t1592 + t150 * t804 + y * t809) + x * (0.5714285715e-1 * t1599 + t187 * t795 + t150 * t809 + y * t821)) + t772 * (t8 * rho_0 * (t1648 + t1799) + rho_0 * t2096)));
        A(4,1) = - (0.7499999998e1 * t783 * t10 * (t125 * t2158 + t772 * t2680));
        A(4,2) = - (0.7499999998e1 * t1554 * (t779 * t772 * (t1417 * t2690 + t1428 * (t2 * t2694 + t2697 - t2699) + t1443 * (t2 * t2708 + t32 * t2704 + t2711 - t2713) + t107 * (t57 * t2704 + t32 * (-0.2000000000e1 * t2698 + 0.2000000000e1 * t2696) + t2 * t2723 + t2726 - t2728) + t72 * (t82 * t2694 + t57 * t2708 + t32 * t2723 + t2 * (0.1333333333e1 * t2725 - 0.1333333333e1 * t2727) + t2739 - t2741) + z * (t377 * t2690 + t82 * (t2697 - t2699) + t57 * (t2711 - t2713) + t32 * (t2726 - t2728) + t2 * (t2739 - t2741) - 0.6666666665e-1 * t148 * t2686 + 0.6666666665e-1 * t137 * t1413)) - 0.3999999999e0 * t1443 * y * mu * t11 + t107 * (-0.1200000000e1 * t1586 + y * (t1496 - t2767)) + 0.7999999998e1 * t181 * t41 * t47 + t72 * (-0.1200000000e1 * t1592 + t150 * (t1506 - 0.2399999999e1 * t2107) + y * (0.9999999998e0 * t2780 - t2782)) + t3 * (-0.1200000000e2 * y * t1 * t41 * t47 - 0.1200000000e2 * t150 * t41 * t47) + z * (-0.3999999999e0 * t1599 + t187 * (-t1527 - t2767) + t150 * (-0.5999999998e1 * t2780 - t2782) + y * (-0.2999999999e1 * t12 * t41 * t40 - 0.3999999999e0 * t2125)) + 0.9999999998e0 * y * t12 * t41 * t47 + 0.2000000000e1 * t150 * t1 * t41 * t47 + 0.9999999998e0 * t187 * t41 * t47) * t10);
        A(4,3) = - t1576;
        A(4,4) = - (t1567 * rho_0 * (t922 - 0.2e1 * t1561 + t122 + t2826 + t1558 + t2827 + t1563));
        A(4,5) = t2838;
        A(5,0) = - (0.7499999998e1 * t782 * t770 * (t125 * (-0.3999999999e0 * z * t813 * t11 + t252 * (-0.1200000000e1 * z * t23 * t11 + t2849 - t2852 - t2854) + t136 * (-0.1200000000e1 * z * t33 * t11 + t2 * t2863 - t2867 - t2869 + t2871) + x * (-0.3999999999e0 * z * t58 * t11 + t32 * t2877 + t2 * t2879 - t2883 + t2885 + t2887)) + t772 * t3252) * t2839);
        A(5,1) = - (0.7499999998e1 * t782 * t770 * (t125 * (-0.3999999999e0 * z * t1598 * t11 + t187 * (-t2852 - t2854 - 0.1200000000e1 * z * t2106 * t11 + t2849) + t150 * (-0.1200000000e1 * z * t2112 * t11 + t1 * t2863 - t2867 - t2869 + t2871) + y * (-0.3999999999e0 * z * t2124 * t11 + t12 * t2877 + t1 * t2879 - t2883 + t2885 + t2887)) + t772 * t3591) * t2839);
        A(5,2) = (t10 * t782 * t770 * (t3653 + t3686 + t3735 + t3775));
        A(5,3) = - (0.5000000000e0 * t1580);
        A(5,4) = t2838;
        A(5,5) = - (t1567 * rho_0 * (t123 + t1560 + t2826 + t1558 - t1561 + t2827 + t1562));



        t1 = (double) x * (double) x;
        t2 = (double) y * (double) y;
        t3 = z * z;
        t5 = sqrt(t1 + t2 + t3);
        double t9 = exp(0.1e1 / H * (-t5 + h_0 + R_e));
        t10 = t9 * rho_0;
        t12 = (double) omega * (double) omega;
        t18 = (double) Vx * (double) Vx;
        t19 = (double) Vy * (double) Vy;
        double t20 = Vz * Vz;
        double t22 = sqrt(t12 * (t1 + t2) + (double) (2 * omega * (Vx * y - Vy * x)) + t18 + t19 + t20);
        double t26 = 0.1e1 / ball / Cd;

        A(3,6) = -0.5e0 * t26 * (double) (omega * y + Vx) * t22 * t10;
        A(4,6) = t26 * t22 * t9 * (0.5e0 * (double) omega * (double) x - 0.5e0 * (double) Vy) * rho_0;
        A(5,6) = -0.5e0 * t26 * Vz * t22 * t10;



        int i;
        int j;

        for (i=0; i < d; i++) {
          for (j=0; j < d; j++) {  
            Phi(i,j) = X[d+d*i+j];
          }
        }
       
        Eigen::MatrixXf dPhidt(d,d);
        dPhidt = A*Phi;
        for (i=0; i < d; i++) {
          for (j=0; j < d; j++) {  
            dXdt[d+d*i+j] = dPhidt(i,j);
          }
        }
      return dXdt;
      }


      state_type SpaceCraft::STM_dot (const state_type X , const double t ){
        // Provides derivative of STM with J2 and J3
        // Input (and output) are (d*d+1) vectors

        state_type dXdt(dimensions*(dimensions+1),0.0);

        double x = X[0];        
        double y = X[1];        
        double z = X[2];        
        int d = dimensions;
 
        double R_e = env.R_e;
        double mu = env.mu;
        double J2;
        double J3;
        if (LEG_GRAV_FIELD_) {
          J2 = - env.Cgrav(2,0);
          J3 = - env.Cgrav(3,0);
        } else { J2 = env.J_2;
          J3 = env.J_3;
        }
        
        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(d,d);
        Eigen::MatrixXf Phi = Eigen::MatrixXf::Zero(d,d);

        double t1 = z * z;
        double t2 = (t1 * t1);
        double t3 = t2 * t2;
        double t5 = (x * x);
        double t7 = (y * y);
        double t9 = R_e * R_e;
        double t10 = (t9 * J2);
        double t11 = 3 * t10;
        double t13 = t1 * t2;
        double t15 = R_e * t9;
        double t16 = (t15 * J3);
        double t17 = (z * t2);
        double t18 = ( t17 * t16);
        double t20 = t5 * t5;
        double t26 = ( t7 * t7);
        double t27 = 3 * t26;
        double t29 = t9 * J2 * t7;
        double t35 = (z * t1);
        double t39 = t5 + t7;
        double t57 = t39 * t39;
        double t66 = t5 + t7 + t1;
        double t67 = t66 * t66;
        double t68 = t67 * t67;
        double t70 = sqrt(t66);
        double t72 = 0.1e1 / t70 / t66 / t68;
        double t78 = 3 * t5;
        double t79 = 3 * t7;
        double t99 = 0.3e1 * mu * t72 * (t13 + (t2 * (-15 * t10 + t78 + t79)) - (35 * t35 * t16) + 0.3e1 * t1 * (-0.25e2 / 0.6e1 * t10 + t5 + t7) * t39 + 0.35e2 / 0.2e1 * z * t39 * t16 + (0.5e1 / 0.2e1 * t10 + t5 + t7) * t57) * x * y;
        double t122 = ( (t35 * t2) + (t17 * (-10 * t10 + t78 + t79)) - (20 * t2 * t16) + 0.3e1 * t35 * (-0.5e1 / 0.6e1 * t10 + t5 + t7) * t39 + 0.30e2 * t1 * t39 * t16 + z * (0.15e2 / 0.2e1 * t10 + t5 + t7) * t57 - 0.5e1 / 0.2e1 * t57 * t16) * mu;
        double t124 = 0.3e1 * t122 * t72 * x;
        double t170 = 0.3e1 * t122 * t72 * y;

        A(0,3) = 1;
        A(1,4) = 1;
        A(2,5) = 1;
        A(3,0) = (0.2e1 * mu * t72 * (-t3 / 0.2e1 + t13 * (- t5 / 0.2e1 - (2 * t7) + t11) + (5 * t18) + t2 * (0.3e1 / 0.2e1 * t20 + t5 * (-0.69e2 / 0.4e1 * t10 - 0.3e1 / 0.2e1 * t7) - t27 + 0.21e2 / 0.4e1 * t29) - 0.205e3 / 0.4e1 * t35 * ( t5 - t7 / 0.41e2) * t16 + 0.5e1 / 0.2e1 * t1 * (t20 + t5 * (-0.69e2 / 0.10e2 * t10 + t7 / 0.5e1) + 0.3e1 / 0.5e1 * t29 - 0.4e1 / 0.5e1 * t26) * t39 + 0.45e2 / 0.2e1 * z * ( t5 - t7 / 0.6e1) * J3 * t39 * t15 + (t20 + t5 * ( t11 + t7 / 0.2e1) - 0.3e1 / 0.4e1 * t29 - t26 / 0.2e1) * t57));
        A(3,1) = t99;
        A(3,2) = t124;
        A(4,0) = t99;
        A(4,1) = - (mu * (t3 + t13 * (-6 * t10 + 4 * t5 + t7) - (10 * t18) + t2 * (0.6e1 * t20 + t5 * (-0.21e2 / 0.2e1 * t10 + t79) - t27 + 0.69e2 / 0.2e1 * t29) - 0.5e1 / 0.2e1 * t35 * (t5 - 41 * t7) * t16 + 0.4e1 * t1 * (t20 + t5 * (-0.3e1 / 0.4e1 * t10 - t7 / 0.4e1) + 0.69e2 / 0.8e1 * t29 - 0.5e1 / 0.4e1 * t26) * t39 + 0.15e2 / 0.2e1 * z * t39 * (t5 - 6 * t7) * t16 + (t20 + t5 * (0.3e1 / 0.2e1 * t10 - t7) - 0.6e1 * t29 - (2 * t26)) * t57) * t72);
        A(4,2) = t170;
        A(5,0) = t124;
        A(5,1) = t170;
        A(5,2) = - (mu * t72 * (-0.2e1 * t3 + t13 * (12 * t10 - 5 * t5 - 5 * t7) + (20 * t18) - 0.3e1 * t2 * (8 * t10 + t5 + t7) * t39 - 0.100e3 * t35 * t39 * t16 + t1 * (-0.63e2 / 0.2e1 * t10 + t5 + t7) * t57 + 0.75e2 / 0.2e1 * z * t57 * t16 + (0.9e1 / 0.2e1 * t10 + t5 + t7) * t39 * t57));

        int i;
        int j;

        for (i=0; i < d; i++) {
          for (j=0; j < d; j++) {  
            Phi(i,j) = X[d+d*i+j];
          }
        }
       
        Eigen::MatrixXf dPhidt(d,d);
        dPhidt = A*Phi;
        for (i=0; i < d; i++) {
          for (j=0; j < d; j++) {  
            dXdt[d+d*i+j] = dPhidt(i,j);
          }
        }
      return dXdt;
      }



      Eigen::Matrix3f SpaceCraft::Rot_axis_(double angle, int ax){

        Eigen::Matrix3f R = Eigen::Matrix3f::Identity(3,3);
        if (ax == 1){
          R(1,1) = std::cos(angle);
          R(1,2) = std::sin(angle);
          R(2,1) = -std::sin(angle);
          R(2,2) = std::cos(angle);
        }
        if (ax == 2){
          R(0,0) = std::cos(angle);
          R(0,2) = -std::sin(angle);
          R(2,0) = std::sin(angle);
          R(2,2) = std::cos(angle);
        }
        if (ax == 3){
          R(0,0) = std::cos(angle);
          R(0,1) = std::sin(angle);
          R(1,0) = -std::sin(angle);
          R(1,1) = std::cos(angle);
        }
      
        return R;

      }


      Eigen::Matrix3f SpaceCraft::orthoDCM_(Eigen::Matrix3f DCM){
        Eigen::Matrix3f eye = Eigen::Matrix3f::Identity(3,3);
        Eigen::Matrix3f delT = DCM*DCM.transpose()- eye;
        double num = 1;
        double den = 2;
        double sign = -1;
        Eigen::Matrix3f PreorthDCM = Eigen::Matrix3f::Identity(3,3);
        Eigen::Matrix3f PreorthDCM_Cumul = Eigen::Matrix3f::Identity(3,3);
        for (int n= 0; n>10; n++) {
          PreorthDCM_Cumul *= sign*num/den*delT;
          PreorthDCM += PreorthDCM_Cumul;
          num += 2.;
          den += 2.;
          sign *= -1;
        }

        Eigen::Matrix3f orthDCM = DCM*PreorthDCM;
        return orthDCM;
      }

      void SpaceCraft::getMisc(double JD_UTC, double* Ome_moon, double* Delta_Psi_1980,
                               double* epsilon_1980, double* epsilon_bar_1980) {
      // get parameters useful for rotation
        double J2000 = 2451545.;
        double MJD_J2000 = J2000 - 2400000.5;
        double MJD_UTC = JD_UTC - 2400000.5;
        double MJD_TAI = MJD_UTC + 37./(3600.*24);
        double MJD_TT = MJD_TAI + 32.184/(3600.*24);
        double T_TT = (MJD_TT - MJD_J2000)/36525.;

        double MJD_finals_start = env.finals_all(0,0);
        int index_finals_all = (int) floor(MJD_UTC - MJD_finals_start);
        double remainder = std::remainder(MJD_UTC, 1.);
        double deg2rad = M_PI/180.;
        double arcsec2rad = M_PI/180./3600.;
        double r = 360.; // 360 
        double T_TT2 = T_TT*T_TT;
        double T_TT3 = T_TT2 * T_TT;
        double DeltaUTC = (1. - remainder)*env.finals_all(index_finals_all,18) + remainder*env.finals_all(index_finals_all+1,18); 

        double epsilon_0 = 23.439291*deg2rad; // rad
        *epsilon_bar_1980 = (23.439291 - 0.0130042*T_TT - 1.64e-7*T_TT2 + 5.04e-7*T_TT3)*deg2rad; // rad

        double M_moon = (134.96298139   + (1325*r  + 198.8673981)*T_TT  + 0.0086972*T_TT2 + 1.78e-5*T_TT3)*deg2rad; // rad
        double M_sun = (357.52772333    + (99*r    + 359.0503400)* T_TT - 0.0001603*T_TT2 - 3.3e-6*T_TT3)*deg2rad; // rad
        double u_M_moon = (93.27191028  + (1342.*r + 82.0175381)*T_TT   - 0.0036825*T_TT2 + 3.1e-6*T_TT3)*deg2rad; // rad
        double D_sun = (297.85036306    + (1236.*r + 307.1114800)*T_TT  - 0.0019142*T_TT2 + 5.3e-6*T_TT3)*deg2rad; // rad
        *Ome_moon = (125.04452222 - (5.*r    + 134.1362608)*T_TT  + 0.0020708*T_TT2 + 2.2e-6*T_TT3)*deg2rad; // rad

        *Delta_Psi_1980 = 0.;
        double Delta_epsilon_1980 = 0.;
        double an1;
        double an2;
        double an3;
        double an4;
        double an5;
        double A;
        double B;
        double C;
        double D;
        double a_P;
      
        for (int i = 0; i < 106; i++) {
          an1 = env.nut80(i,0);
          an2 = env.nut80(i,1);
          an3 = env.nut80(i,2);
          an4 = env.nut80(i,3);
          an5 = env.nut80(i,4);
          A   = env.nut80(i,5);
          B   = env.nut80(i,6);
          C   = env.nut80(i,7);
          D   = env.nut80(i,8);

          a_P = an1 * M_moon + an2*M_sun + an3*u_M_moon + an4*D_sun + an5*(*Ome_moon);

          *Delta_Psi_1980 += (A + B * T_TT) * std::sin(a_P);
          Delta_epsilon_1980 += (C + D * T_TT) * std::cos(a_P);
        }

        *Delta_Psi_1980 *= 1e-4*arcsec2rad;
        Delta_epsilon_1980 *= 1e-4*arcsec2rad;
        *epsilon_1980 = *epsilon_bar_1980 + Delta_epsilon_1980;

      }

      Eigen::Matrix3f SpaceCraft::getNutationMat(double JD_UTC) {


        double Ome_moon;
        double Delta_Psi_1980;
        double epsilon_1980;
        double epsilon_bar_1980;
        getMisc(JD_UTC, &Ome_moon, &Delta_Psi_1980, &epsilon_1980, &epsilon_bar_1980);

        Eigen::Matrix3f N = Rot_axis_(- epsilon_bar_1980,1)*Rot_axis_(Delta_Psi_1980,3)*Rot_axis_(epsilon_1980,1);
        return N;

      }


      Eigen::Matrix3f SpaceCraft::getPrecessionMat(double JD_UTC){

        double J2000 = 2451545.;
        double MJD_J2000 = J2000 - 2400000.5;
        double MJD_UTC = JD_UTC - 2400000.5;
        double MJD_TAI = MJD_UTC + 37./(3600.*24);
        double MJD_TT = MJD_TAI + 32.184/(3600.*24);
        double T_TT = (MJD_TT - MJD_J2000)/36525.;

        double T_TT2 = T_TT*T_TT;
        double T_TT3 = T_TT2 * T_TT;

        double arcsec2rad = M_PI/180./3600.;

        double xi_A =    (2306.2181*T_TT + 0.30188*T_TT2   + 0.017998*T_TT3)*arcsec2rad;
        double theta_A = (2004.3109*T_TT - 0.426625*T_TT2 - 0.041833*T_TT3)*arcsec2rad;
        double z_A =     (2306.2181*T_TT + 1.09468*T_TT2   + 0.018203*T_TT3)*arcsec2rad;
        
        Eigen::Matrix3f P = Rot_axis_(xi_A,3)*Rot_axis_(-theta_A,2)*Rot_axis_(z_A,3);
        return P;
      }


      double SpaceCraft::getAlpha_G(double JD_UTC) {
        double arcsec2rad = M_PI/180./3600.;
        double Ome_moon;
        double Delta_Psi_1980;
        double epsilon_1980;
        double epsilon_bar_1980;
        getMisc(JD_UTC, &Ome_moon, &Delta_Psi_1980, &epsilon_1980, &epsilon_bar_1980);

        double J2000 = 2451545.;
        double MJD_J2000 = J2000 - 2400000.5;
        double MJD_UTC = JD_UTC - 2400000.5;

        double MJD_finals_start = env.finals_all(0,0);
        int index_finals_all = (int) floor(MJD_UTC - MJD_finals_start);
        double remainder = std::remainder(MJD_UTC, 1.);

        double DeltaUTC = (1. - remainder)*env.finals_all(index_finals_all,18) + remainder*env.finals_all(index_finals_all+1,18); 
        double MJD_UT1 = MJD_UTC + DeltaUTC/(24.*3600);
        double DeltaT = MJD_UT1 - MJD_J2000;

        double Theta_GMST_UT1 = 4.894961212823058751375704430 + (6.3003880989848935522776513720 + (5.07520999411359147805523e-15 - 9.253097568194335640067190688e-24*DeltaT)*DeltaT)*DeltaT;
        double Eq_equinox =  Delta_Psi_1980*std::cos(epsilon_bar_1980) + (0.00264*std::sin(Ome_moon) + 0.000063*std::sin(2*Ome_moon))*arcsec2rad;
        double alpha_G = Theta_GMST_UT1 + Eq_equinox; 
        return alpha_G;
      }

      Eigen::Matrix3f SpaceCraft::getPolarWMat(double JD_UTC) {
        double MJD_UTC = JD_UTC - 2400000.5;
        double MJD_finals_start = env.finals_all(0,0);
        int index_finals_all = (int) floor(MJD_UTC - MJD_finals_start);
        double remainder = std::remainder(MJD_UTC, 1.);
        double arcsec2rad = M_PI/180./3600.;

        double x_p = ((1 - remainder)*env.finals_all(index_finals_all,16) + remainder*env.finals_all(index_finals_all+1,16))*arcsec2rad;
        double y_p = ((1 - remainder)*env.finals_all(index_finals_all,17) + remainder*env.finals_all(index_finals_all+1,17))*arcsec2rad;
        Eigen::Matrix3f W = Rot_axis_(y_p,1) * Rot_axis_(x_p,2);
        return W;
      }

      Eigen::Matrix3f SpaceCraft::getS_primeMat(double JD_UTC) {

        double alpha_G = getAlpha_G(JD_UTC);
        Eigen::Matrix3f S_prime = Rot_axis_(-alpha_G, 3);
        return S_prime;
      }


      Eigen::Matrix3f SpaceCraft::getS_prime_dotMat(double JD_UTC) {

        double alpha_G = getAlpha_G(JD_UTC);
        Eigen::Matrix3f S_prime_dot;
        double sa = std::sin(alpha_G);
        double ca = std::cos(alpha_G);
        S_prime_dot << -sa, ca, 0.0,
                       -ca, -sa, 0.0,
                       0.0, 0.0, 0.0;
        return env.omega * S_prime_dot;
      }

      Eigen::Matrix3f SpaceCraft::getECEF2ECI_dot(double JD_UTC){

        Eigen::Matrix3f P = getPrecessionMat(JD_UTC);
        Eigen::Matrix3f N = getNutationMat(JD_UTC);
        Eigen::Matrix3f S_prime_dot = getS_prime_dotMat(JD_UTC);
        Eigen::Matrix3f W = getPolarWMat(JD_UTC);
        Eigen::Matrix3f T_dot = P*N*S_prime_dot*W;
       
        return orthoDCM_(T_dot);
      }

      Eigen::Matrix3f SpaceCraft::getECEF2ECI(double JD_UTC){
        // generate rotation matrix from ECI to ECEF
        // T = WSNP  

        // Time 
        Eigen::Matrix3f P = getPrecessionMat(JD_UTC);
        Eigen::Matrix3f N = getNutationMat(JD_UTC);
        Eigen::Matrix3f S_prime = getS_primeMat(JD_UTC);
        Eigen::Matrix3f W = getPolarWMat(JD_UTC);
        Eigen::Matrix3f T = P*N*S_prime*W;
       
        return orthoDCM_(T);
        
      } 

      void SpaceCraft::propagate(double dt){
        if (LEG_GRAV_FIELD_) {
          env.J_2 = 0.;
        }
        int d = dimensions;
        state_type x(d*(d+1), 0.0);
        int i = 0;
        for (auto & pos: cart.ipos) {
          x[i] = cart.ipos[i];
          x[3+i] = cart.ivel[i];
          i++;
        }
        for (i=0; i < d; i++) {
          for (int j=0; j < d; j++) {  
            x[d+d*i+j] = cart.STM(i,j);
          }
        }
        double init_dt;
        if (std::abs(dt) < std::abs(integ.dtGuess)) {
          init_dt = dt;
        } else { init_dt = integ.dtGuess;
        }
        if (dt < 0) {
          init_dt = - std::abs(init_dt);
        } else { init_dt = std::abs(init_dt);
        }
        ode::integrate_adaptive( ode::make_controlled( integ.absTol , integ.relTol ,
			error_stepper_type() ) , *this , x , t , t+dt , init_dt );
        i = 0;
        for (auto & pos: cart.ipos) {
          cart.ipos[i] = x[i];
          cart.ivel[i] = x[i+3];
          i++;
        }
        for (i=0; i < d; i++) {
          for (int j=0; j < d; j++) {  
            cart.STM(i,j) = x[d+d*i+j];
          }
        }



        t += dt;
        iCart2Kepl_();
  
      }

      void SpaceCraft::setGravField (int degree, int order, std::string Cgrav_name, std::string Sgrav_name, double JD_UTC_0, bool normalized) {
        // loads gravity field of degree and order
        LEG_GRAV_FIELD_ = true;
        env.degree = degree;
        env.order = order;
        env.JD_UTC_0 = JD_UTC_0;
        env.Cgrav = LoadDatFile(Cgrav_name, degree+1, env.order+1);
        env.Sgrav = LoadDatFile(Sgrav_name, degree+1, env.order+1);
        // Nutation and rotation data is hardcoded
        env.nut80 = LoadDatFile("nut80.csv", 106, 10);
        env.finals_all = LoadDatFile("finals_all_mod_red.csv", 33, 21);
        
        double normalization;
        if (normalized) {
          for (int n = 0; n < env.degree+1; n++){ // unnormalize and introduce (-1)^m consistent with boost 
            normalization = std::sqrt(((2.0*n + 1.0)));
            env.Cgrav(n,0) *= normalization;
            env.Sgrav(n,0) *= normalization;
            for (int m = 1; m < std::min(env.order+1,n+1); m++){
              normalization = std::sqrt((2.0*(2.0*n + 1.0)
                *boost::math::factorial<double>(n - m)/boost::math::factorial<double>(n + m)));
              env.Cgrav(n,m) *= normalization;
              env.Sgrav(n,m) *= normalization;
            }
          }
        } else { // only introduce (-1)^m consistent with boost 
          for (int n = 0; n < env.degree+1; n++){ // unnormalize and introduce (-1)^m consistent with boost 
            for (int m = 0; m < std::min(env.order+1, n+1); m++){
              env.Cgrav(n,m) /= std::pow(-1.,m);
              env.Sgrav(n,m) /= std::pow(-1.,m);
            }
          }
        }
      }

      void SpaceCraft::setLuniSolar (double AU, double R_e, double mu_Sun, double mu_Moon, double JD_UTC_0) {
        LUNISOLAR_ = true;
        env.AU = AU;
        env.R_e = R_e;
        env.mu_Sun = mu_Sun;
        env.mu_Moon = mu_Moon;
        env.JD_UTC_0 = JD_UTC_0;

      }

      Eigen::Vector3f SpaceCraft::getLuniSolarAcc(const state_type &X, const double t){
        double JD_UTC = env.JD_UTC_0 + t/86400.0;
        Eigen::Vector3f ipos;
        ipos(0) = X[0];
        ipos(1) = X[1];
        ipos(2) = X[2];
        Eigen::Vector3f pos_Moon = SpaceCraft::MoonPosition(JD_UTC);
        Eigen::Vector3f pos_Sun = SpaceCraft::SunPosition(JD_UTC, JD_UTC);
        Eigen::Vector3f r_satMoon = pos_Moon - ipos; 
        Eigen::Vector3f r_satSun = pos_Sun - ipos; 
        
        Eigen::Vector3f lsacc = Eigen::Vector3f::Zero(3);
        double r_sMn = r_satMoon.norm();
        double r_Mn = pos_Moon.norm();
        double r_sSn = r_satSun.norm();
        double r_Sn = pos_Sun.norm();
        lsacc += env.mu_Moon*(r_satMoon/(r_sMn*r_sMn*r_sMn)
          - pos_Moon/(r_Mn*r_Mn*r_Mn));
        lsacc += env.mu_Sun*(r_satSun/(r_sSn*r_sSn*r_sSn)
          - pos_Sun/(r_Sn*r_Sn*r_Sn));
       
        return lsacc;
      }
      void SpaceCraft::setStations ( std::vector<Eigen::Vector3f> stations_pos_ECEF,
                                     double JD_UTC_0, double omega) {
        env.nut80 = LoadDatFile("nut80.csv", 106, 10);
        env.finals_all = LoadDatFile("finals_all_mod_red.csv", 33, 21);
        env.trueMeas = LoadDatFile("measurements6days.csv", 2571,4);
        env.JD_UTC_0 = JD_UTC_0;
        env.stations_pos_ECEF = stations_pos_ECEF;
        env.omega = omega;
      }

      double SpaceCraft::LightTimeCorrection(Eigen::Vector3f station_pos_ECI,
                                           Eigen::Vector3f* corrected_sat_pos,
                                           Eigen::Vector3f* corrected_sat_vel) {
        double lt = (cart.ipos - station_pos_ECI).norm()/env.c; 
        double sat_r = cart.ipos.norm();
        Eigen::Vector3f approx_accel = - env.mu * cart.ipos / (sat_r*sat_r*sat_r);  
        *corrected_sat_pos = cart.ipos - cart.ivel * lt - 0.5 * approx_accel*lt*lt;
        *corrected_sat_vel = cart.ivel - approx_accel * lt;
        return lt;
        
      }

      Eigen::MatrixXf SpaceCraft::getHmatrix(int station, double t) {
 
        Eigen::MatrixXf H = Eigen::MatrixXf::Zero(2,7);

        double JD_UTC = env.JD_UTC_0 + t/86400.;

        Eigen::Vector3f pos = cart.ipos;
        Eigen::Vector3f vel = cart.ivel;
        Eigen::Matrix3f T = getECEF2ECI(JD_UTC);
        Eigen::Vector3f stat_pos = T * env.stations_pos_ECEF[station];


        
        Eigen::Matrix3f P = getPrecessionMat(JD_UTC);
        Eigen::Matrix3f N = getNutationMat(JD_UTC);
        Eigen::Matrix3f S_prime = getS_primeMat(JD_UTC);
        Eigen::Matrix3f W = getPolarWMat(JD_UTC);
        
        Eigen::Vector3f r_W = orthoDCM_(W)*env.stations_pos_ECEF[station];
        Eigen::Vector3f omegaVec = Eigen::Vector3f::Zero(3);
        omegaVec(2) = env.omega; 

        Eigen::Vector3f stat_vel = orthoDCM_(P*N*S_prime)* (omegaVec.cross(r_W)); 
                 

        double t1 = pos(0) - stat_pos(0);
        double t2 = ( pos(0) * pos(0));
        double t5 = ( pos(1) * pos(1));
        double t8 = ( pos(2) * pos(2));
        double t11 = ( stat_pos(0) * stat_pos(0));
        double t12 = ( stat_pos(1) * stat_pos(1));
        double t13 = ( stat_pos(2) * stat_pos(2));
        double t14 = -2 * pos(0) * stat_pos(0) - 2 * pos(1) * stat_pos(1) - 2 * pos(2) * stat_pos(2) + t11 + t12 + t13 + t2 + t5 + t8;
        double t15 = sqrt( t14);
        double t16 = 0.1e1 / t15;
        double t17 = t16 * t1;
        double t19 = t16 * (pos(1) - stat_pos(1));
        double t21 = t16 * (pos(2) - stat_pos(2));
        double t22 = (vel(0) - stat_vel(0));
        double t25 = stat_vel(1) - vel(1);
        double t26 = t1 * t25;
        double t31 = stat_vel(2) - vel(2);
        double t32 = t1 * t31;
        double t37 = -stat_pos(2) * t22;
        double t42 = 0.1e1 / t15 / t14;
        double t47 = -t22 * pos(1);
        double t51 = pos(0) * t22;
        double t56 = - pos(1) * t31;
        double t57 = t25 * pos(2);
        double t63 = t56 + t57;
        double t69 = -t22 * pos(2);
        H(0,0) = t17;
        H(0,1) = t19;
        H(0,2) = t21;
        H(1,0) = t42 * ( (t5 * t22) + pos(1) * (- (2 * stat_pos(1) * t22) + t26) + (t8 * t22) + pos(2) * (- (2 * stat_pos(2) * t22) + t32) + (t12 * t22) - stat_pos(1) * t26 - stat_pos(2) * ( t37 + t32));
        H(1,1) = t42 * (- t11 * t25 + stat_pos(0) * (- (stat_pos(1) * t22) + 0.2e1 * pos(0) * t25 - t47) + stat_pos(1) * (-t31 * pos(2) + stat_pos(2) * t31 + t51) - t13 * t25 + stat_pos(2) * (t56 + 0.2e1 * t57) - t2 * t25 + (pos(0) * t47) - pos(2) * t63);
        H(1,2) = t42 * (- t11 * t31 + stat_pos(0) * (0.2e1 * pos(0) * t31 + t37 - t69) - t12 * t31 + stat_pos(1) * (t25 * stat_pos(2) + 0.2e1 * pos(1) * t31 - t57) + stat_pos(2) * (-t25 * pos(1) + t51) - t2 * t31 + (pos(0) * t69) + pos(1) * t63);
        H(1,3) = t17;
        H(1,4) = t19;
        H(1,5) = t21;
     
      return H;
      }

      std::vector<double> SpaceCraft::cart2meas(Eigen::Vector3f corrected_sat_pos,
                                     Eigen::Vector3f corrected_sat_vel,
                                     Eigen::Vector3f station_pos_ECI,
                                     Eigen::Vector3f station_vel_ECI) {
        // index 0: range
        // index 1: range rate
        // index 2: right ascension
        // index 3: declination
        // index 4: right ascension rate
        // index 5: declination rate
        std::vector<double> expmeas;
        Eigen::Vector3f rho_vec = corrected_sat_pos - station_pos_ECI; 
        Eigen::Vector3f rho_dot_vec = corrected_sat_vel - station_vel_ECI; 
        double rho_hor = std::sqrt(rho_vec(0)*rho_vec(0) + rho_vec(1)*rho_vec(1));
        expmeas.push_back(rho_vec.norm()); // range
        expmeas.push_back((rho_vec.dot(rho_dot_vec))/expmeas[0]); // range_rate
        if (std::abs(rho_hor) > 1e-12) {
          expmeas.push_back(std::atan2(rho_vec(1), rho_vec(0)));
        } else { expmeas.push_back(std::atan2(rho_dot_vec(1), rho_dot_vec(0)));
        }
        expmeas.push_back(std::asin(rho_vec(2)/expmeas[0]));
        expmeas.push_back(rho_dot_vec(0)*rho_vec(1) - rho_dot_vec(1)*rho_vec(0)/(-rho_hor*rho_hor));
        expmeas.push_back((rho_dot_vec(2) - expmeas[1]*std::sin(expmeas[3]))/rho_hor);
        return expmeas;  
      }

      double SpaceCraft::getExpectedMeasurement(int station, std::vector<int> indices,
                                        double t, std::vector<double>* expmeas) {
        // outputs the light time delay
        // index 0: range
        // index 1: range rate
        // index 2: right ascension
        // index 3: declination
        // index 4: right ascension rate
        // index 5: declination rate
        double JD_UTC = env.JD_UTC_0 + t/86400;

        Eigen::Matrix3f T = getECEF2ECI(JD_UTC);
        Eigen::Vector3f station_pos_ECI = T * env.stations_pos_ECEF[station];

        Eigen::Vector3f corrected_sat_pos;
        Eigen::Vector3f corrected_sat_vel;
        double lt_corr = LightTimeCorrection(station_pos_ECI, &corrected_sat_pos, &corrected_sat_vel );

        Eigen::Matrix3f P = getPrecessionMat(JD_UTC);
        Eigen::Matrix3f N = getNutationMat(JD_UTC);
        Eigen::Matrix3f S_prime = getS_primeMat(JD_UTC);
        Eigen::Matrix3f W = getPolarWMat(JD_UTC);
        Eigen::Vector3f r_W = orthoDCM_(W)*env.stations_pos_ECEF[station];
        Eigen::Vector3f omegaVec = Eigen::Vector3f::Zero(3);
        omegaVec(2) = env.omega; 

        Eigen::Vector3f station_vel_ECI = orthoDCM_(P*N*S_prime)* (omegaVec.cross(r_W)); 
        

        std::vector<double> radar = cart2meas(corrected_sat_pos, corrected_sat_vel,
                          station_pos_ECI, station_vel_ECI); 
        bool bias = true;
        if (bias) {
          if (station == 2) {
            radar[0] += 20.571e-3; 
          }
        }
        for (int n=0; n < indices.size(); n++){
          (*expmeas).push_back(radar[indices[n]]);
        }

        return lt_corr;
      }

      void SpaceCraft::setSRPCannon ( double solarRad, double srp_ball, double AU,
                                      double R_e, double JD_UTC_0, double alpha_pen,
                                      double alpha_umb) {

        SRPCANNON_ = true;
        env.p_srp = solarRad/(env.c*1e3)/1e3; // kN/m^2
        env.AU = AU;
        env.R_e = R_e;
        env.alpha_pen = alpha_pen;
        env.alpha_umb = alpha_umb;
        env.JD_UTC_0 = JD_UTC_0;
        env.srp_ball = srp_ball; // "ballistic" coefficient for SRP (C_R * A / m)
      }

      Eigen::Vector3f SpaceCraft::getSRPCannon (const state_type &X, const double t) {
        double JD_UTC = env.JD_UTC_0 + t/86400.;
        Eigen::Vector3f ipos;
        ipos(0) = X[0];
        ipos(1) = X[1];
        ipos(2) = X[2];
        Eigen::Vector3f pos_sun = SpaceCraft::SunPosition(JD_UTC, JD_UTC);
        Eigen::Vector3f r_satSun = pos_sun - ipos; 

        double shadow = ratioSun(pos_sun, ipos); 
        Eigen::Vector3f SRP_accel = - shadow * env.p_srp/env.srp_ball*r_satSun/r_satSun.norm();

        return SRP_accel;
      }

      double SpaceCraft::ratioSun(Eigen::Vector3f pos_Sun, Eigen::Vector3f pos_sat) {
        // Return 0 if in the shadow, 1 if in light, in-between if penumbra
        // Algo 34 page 301 from Vallado
        double ratioSun = 1.0;
        double dotprod = pos_Sun.transpose()*pos_sat; 
        if (dotprod < 0.0) {
          if (pos_sat.norm() > 60e3) {
            std::cout << "shadow function breaks for very high altitude satellites" << std::endl;
            exit(0);
          }
          double angle = -dotprod/pos_Sun.norm()/pos_sat.norm();
          double sat_horiz = pos_sat.norm() * std::cos(angle);
          double sat_vert = pos_sat.norm() * std::sin(angle);
          double x = env.R_e / std::sin(env.alpha_pen); 
          double pen_vert = std::tan(env.alpha_pen)*(x+sat_horiz); 
          if (sat_vert < pen_vert) {
            double y = env.R_e/std::sin(env.alpha_umb);
            double umb_vert = std::tan(env.alpha_umb)*(y-sat_horiz);
            if (sat_vert < umb_vert) {
              ratioSun = 0.0;
            } else {
              ratioSun = (sat_vert - umb_vert)/(pen_vert - umb_vert);
              //std::cout << ratioSun << std::endl;
              if (ratioSun > 1) {
                exit(0);
              }
              if (ratioSun < 0) {
                exit(0);
              }
            }
          }
          
        }
        return ratioSun;

      }

      void SpaceCraft::setExpDrag (double rho_0, double r_0, double H, double omega, double mass, double Area, double Cd) {

        EXP_DRAGCANNON_ = true;

        env.r_0 = r_0;
        env.rho_0 = rho_0;
        env.H = H;
        env.mass = mass;
        env.Area = Area;
        env.Cd = Cd;
        env.omega = omega;
      
      }
 
      void SpaceCraft::operator() (const state_type &X , state_type &dXdt , const double /* t */ ){
        // This function provides the derivatives
        int d = dimensions;

        dXdt[0] = X[3];
        dXdt[1] = X[4];
        dXdt[2] = X[5];
        Eigen::Vector3f ipos;
        double x = X[0];        
        double y = X[1];        
        double z = X[2];     
        ipos(0) = x;
        ipos(1) = y;
        ipos(2) = z;   

        double t2 = env.R_e * env.R_e;
        double t3 = (t2 * env.J_2);
        double t4 = (x * x);
        double t7 = (y * y);
        double t10 = (z * z);
        double t13 = (t4 * t4);
        double t19 = (t7 * t7);
        double t23 = (t10 *t10);
        double t26 = t4 + t7 + t10;
        double t27 = t26 * t26;
        double radius = std::sqrt(t26);
        double t31 = 1. / radius / t26 / t27;
        double t32 = t31 * (t10*(-12. * t3 + 4. * t4 + 4. * t7) + 3. * t4 * t3 + 3. * t7 * t3 + 4. * t7 * t4 + 2. * t13 + 2. * t19 + 2. * t23);
        dXdt[3] = -t32 * env.mu * x / 0.2e1;
        dXdt[4] = -t32 * env.mu * y / 0.2e1;
        dXdt[5] = -4.5 * z * env.mu * t31 * (0.2e1 / 0.9e1 * t23 + t10 * (-0.2e1 / 0.3e1 * t3 + 0.4e1 / 0.9e1 * t4 + 0.4e1 / 0.9e1 * t7) + (t3 + 0.2e1 / 0.9e1 * t4 + 0.2e1 / 0.9e1 * t7) * (t4 + t7));


        if (EXP_DRAGCANNON_) {
          double rho = env.rho_0 * std::exp(-(radius - env.r_0)/env.H); 
          Eigen::Vector3f V_A;
          Eigen::Vector3f accel_drag;
          V_A[0] = X[3] + env.omega*X[1]; 
          V_A[1] = X[4] - env.omega*X[0]; 
          V_A[2] = X[5];
          accel_drag = -0.5 * rho* V_A * V_A.norm() / env.mass * env.Cd *env.Area;
          dXdt[3] += accel_drag[0];
          dXdt[4] += accel_drag[1];
          dXdt[5] += accel_drag[2];

        }
        if (LEG_GRAV_FIELD_) {
          Eigen::Vector3f dUdspher = Eigen::Vector3f::Zero(3); // [in Zenith East Norht]
                                     // still Cartesian though, but in direction of angles and radius
         
          Eigen::Matrix3f T = getECEF2ECI(env.JD_UTC_0 + t/86400.);
          Eigen::Vector3f rpos = T.inverse()*ipos;
          //std::cout << "ipos" << std::endl;
          //std::cout << ipos << std::endl;
          //std::cout << "rpos" << std::endl;
          //std::cout << rpos << std::endl;
          Eigen::Vector3f rposspher = Cart2spher_vec_(rpos);
          double r = rposspher(0);
          double lambda = rposspher(1);
          double phi = rposspher(2);
          //std::cout << "spher_pos" << std::endl;
          //std::cout << rposspher << std::endl;
          
          double cosmlambda; // sin(mlambda)
          double sinmlambda; // sin(mlambda)
          double sinphi = std::sin(phi); // sin(phi)
          double cosphi = std::cos(phi); // cos(phi)
          double secphi = 1./cosphi; // sec(phi)
          double legnm; // associated legendre poly
          double legnm_1; // associated legendre poly for order n-1
          double cosdlegdphi; // cos(phi) d(P(sphi))/d(sphi)
          double Rrn = (env.R_e/r); // (R/r)^l
          double Rr = (env.R_e/r); // (R/r)^l
          double mur = (env.mu/r); // mu/R
          double mur2 = mur/r; // mu/R^2 

          for (int n=2; n < env.degree + 1; n++){
            Rrn *= Rr;
            for (int m=0; m < std::min(env.order+1,n+1); m++){
              cosmlambda = std::cos(m*lambda);
              sinmlambda = std::sin(m*lambda);
              legnm = std::pow(-1.0,m)*boost::math::legendre_p(n,m,sinphi);
              //std::cout << "legnm" << std::endl;
              //std::cout << legnm << std::endl;
              //std::cout << "s and C" << std::endl;
              //std::cout << env.Sgrav(n,m) << std::endl;
              //std::cout << env.Cgrav(n,m) << std::endl;
              dUdspher(0) += Rrn*(n+1)*legnm*(env.Cgrav(n,m)*cosmlambda + env.Sgrav(n,m)*sinmlambda);
              //std::cout << "dUdspher_0" << std::endl;
              //std::cout << dUdspher(0) << std::endl;

              dUdspher(1) += Rrn*m*legnm*(env.Sgrav(n,m)*cosmlambda - env.Cgrav(n,m)*sinmlambda);
              if (m==n) {
                legnm_1 = 0.;
              } else {
                //legnm_1 = std::pow(-1.0,m)* boost::math::legendre_p(n-1,m,sinphi);
                legnm_1 = std::pow(-1.0,m+1)* boost::math::legendre_p(n,m+1,sinphi);
              }
              //cosdlegdphi = -n*sinphi*secphi*legnm + (n+m)*secphi*legnm_1;
              cosdlegdphi = legnm_1 - m*std::tan(phi)*legnm;
              dUdspher(2) += Rrn*cosdlegdphi*(env.Cgrav(n,m)*cosmlambda + env.Sgrav(n,m)*sinmlambda);

            }
          }
          dUdspher(0) *= (-mur2);
          dUdspher(1) *= (mur)/r/cosphi;
          dUdspher(2) *= (mur)/r;
          //std::cout << "dUdspher" << std::endl;
          //std::cout << dUdspher << std::endl;
     
          Eigen::Vector3f rcart_acc = Rot_axis_(-lambda, 3)*Rot_axis_(phi,2)*dUdspher;
          //std::cout << "rcart_acc" << std::endl;
          //std::cout << rcart_acc << std::endl;
          Eigen::Vector3f icart_acc = T*rcart_acc;
          
          //std::cout << "icart_acc" << std::endl;
          //std::cout << icart_acc << std::endl;

          dXdt[3] += icart_acc[0];
          dXdt[4] += icart_acc[1];
          dXdt[5] += icart_acc[2];
          //std::cout << " " << std::endl;
          //std::cout << -0.0075686041124269562 - dXdt[3] << std::endl;
          //std::cout <<-0.0017471361235711653 -  dXdt[4] << std::endl;
          //std::cout <<-1.4317417903930242e-05 - dXdt[5] << std::endl;
          //exit(0);
        }

        if (LUNISOLAR_) {
          Eigen::Vector3f lsacc = getLuniSolarAcc(X, t);
          dXdt[3] += lsacc(0);
          dXdt[4] += lsacc(1);
          dXdt[5] += lsacc(2);
        }

        if (SRPCANNON_) {
          Eigen::Vector3f SRPacc = getSRPCannon(X, t);
          dXdt[3] += SRPacc(0);
          dXdt[4] += SRPacc(1);
          dXdt[5] += SRPacc(2);
        }
        if (STM_PROP_) {
          state_type dSTMdt = STM_dot2(X, t);
          for (int jj=d; jj < (d*(d+1)); jj++) {
            dXdt[jj] = dSTMdt[jj];
          } 
        }
      }

} 
