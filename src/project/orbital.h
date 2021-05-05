#ifndef OD_SRC_PROJECT_ORBITAL_H_
#define OD_SRC_PROJECT_ORBITAL_H_

#include<iostream>
#include<cmath>
#include<boost/numeric/odeint.hpp>
#include<fstream>
#include<string>
// #include"../../../../eigen/Eigen/Eigen"
#include <Eigen/Dense> //matricies and vectors
#include<boost/array.hpp>

namespace orbital{

  namespace ode = boost::numeric::odeint;
  typedef std::vector< double > state_type;
  typedef ode::runge_kutta_dopri5< state_type > error_stepper_type;
  typedef ode::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;



  // Spacecraft state, Cartesian
  struct Cart{
    Eigen::Vector3f ipos; // inertial position, wrt ECI
    Eigen::Vector3f ivel; // inertial velocity, wrt ECI
    Eigen::Vector3f rpos; // position wrt ECEF
    Eigen::Vector3f rvel; //velocity wrt ECEF
    Eigen::MatrixXf STM; // State Transition Matrix
    Eigen::MatrixXf Cov; // Uncertainty Covariance Matrix
  };


  // Spacecraft state, Keplerian osculating
  struct OE{
    double sma; // osculating orbit semi-major axis
    double ecc; // osculating orbit eccentricity
    Eigen::Vector3f ee; // spacecraft eccentricity vector
    double incl; // osculating orbit inclination
    double RAAN; // osculating orbit right ascension of ascending node
    double ome; // osculating orbit argument of periapsis
    double theta; // osculating orbit true anomaly
    double E; // osculating orbit eccentric anomaly
    double M; // osculating orbit mean anomaly
    double n; // osculating orbit mean angular velocity
    double T_p; // osculating orbit time of perigee passage
  };

  // Environment (and spacecraft)
  struct Env{
    double mu; // gravitational parameter of central body [km^3/s^2]
    double R_e; // planet Equatorial radius [km]
    double J_2; // J_2 [-]
    double J_3; // J_2 [-]
    double AU; // Astronomical Unit [km] (used in algo for Sun position)
    double mu_Moon; // Gravitational parameter of the moon
    double mu_Sun; // Gravitational parameter of the sun

    Eigen::MatrixXf Sgrav; // S coefficients of gravity field
    Eigen::MatrixXf Cgrav; // S coefficients of gravity field
    Eigen::MatrixXf finals_all; // finals (for rotations)
    Eigen::MatrixXf nut80; // nutation data
    int degree; // degree of gravity field
    int order; // order of gravity field
    double JD_UTC_0; // propagation start time

    double omega; // planet rotation rate  [rad/s]
    double rho_0; // density at r_0 [kg/m^2/km]
    double r_0; // altitude at which dnesity is rho_0 [km]
    double H; // density scale height [km]

    double drag_ball; // drag ballistic coefficient of the spacecraft [m/kg]
    double Cd;
    double frontalArea;
    double solarPanelArea;
    double srp_ball; // srp "ballistic" coefficient of the spacecraft [m/kg]
    double alpha_pen; // angle penumbra
    double alpha_umb; // angle umbra
    double mass; // spacecraft mass (useful if using not cannonball models) [kg]

 
    double c = 299792.458; // lightspeed, km/s
    double p_srp; // solar radiation pressure of the Sun

    // Measurements
    std::vector<Eigen::Vector3f> stations_pos_ECEF;
    Eigen::MatrixXf trueMeas; 
  };

  // Variables for integration
  struct Integrator{
    double dtGuess; // initial guess for step size
    double absTol;
    double relTol;
  };

  class SpaceCraft{
    public:

      // Spacecraft state, Cartesian
      Cart cart;
      // Spacecraft state, Keplerian osculating
      OE oe;

      double sigma_Q; // Uncertainty in force model (km/s^2) 

      // Others
      Eigen::Vector3f hh; // spacecraft angular momentum vector
      double h; // spacecraft angular momentum
      double energy_pm; // spacecraft energy (potential considering only point mass)
      double energy_J2; // spacecraft energy (potential considering only point mass and J2)
      double period; // orbital period

      // Time
      double t; // time
      double last_update_time;
      
      // Environment
      Env env;

      // 
      int dimensions;

      // Variables for integration
      Integrator integ;
    
      // Constructor
      SpaceCraft(int d);

      Eigen::VectorXf getKepl();

      Eigen::VectorXf getCartInert();

      Eigen::VectorXf getCartRot();
 
      void setKepl(Eigen::VectorXf Kepl);

      void setCartInert(Eigen::Vector3f ipos, Eigen::Vector3f ivel);

      void setCartInertCov(Eigen::MatrixXf cov);

      void setCartRot(Eigen::Vector3f rpos, Eigen::Vector3f rvel);

      void setExpDrag(double rho_0, double r_0, double H, double omega, double mass, double frontalArea,  double solarPanelArea, double Cd);

      void setGravField (int degree, int order, std::string Cgrav_name, std::string Sgrav_name, double JD_UTC_0, bool normalized);

      void setLuniSolar (double AU, double R_e, double mu_Sun, double mu_Moon, double JD_UTC_0);

      state_type STM_dot (const state_type X , const double t );

      state_type STM_dot2 (const state_type X , const double t );

      Eigen::MatrixXf getHmatrix(int station, double t);

      void KalmanUpdate(std::vector<double> measurement, int station);

      void setSRPCannon (double solarRad, double srp_ball, double AU, double R_e, double JD_UTC_0, 
                         double alpha_pen, double alpha_umb);

      double ratioSun(Eigen::Vector3f pos_Sun, Eigen::Vector3f pos_sat); 

      Eigen::Vector3f getSRPCannon (const state_type &X, const double t);

      void setStations (std::vector<Eigen::Vector3f> stations_pos_ECEF, double JD_UTC_0, double omega);

      void setSTMprop();

      double LightTimeCorrection(Eigen::Vector3f station_pos_ECI,
                               Eigen::Vector3f* corrected_sat_pos,
                               Eigen::Vector3f* corrected_sat_vel);

      double getExpectedMeasurement(int station, std::vector<int> indices,
                            double t, std::vector<double>* expmeas); 

      std::vector<double> cart2meas(Eigen::Vector3f corrected_sat_pos,
                                     Eigen::Vector3f corrected_sat_vel,
                                     Eigen::Vector3f station_pos_ECI,
                                     Eigen::Vector3f station_vel_ECI);

      Eigen::Vector3f getLuniSolarAcc(const state_type &X, const double t);

      // () operator function reserved for orbit propagation with odeint
      void operator()( const state_type &X , state_type &dXdt , const double /* t */ );
   
      void propagate(double dt);

      void getMisc(double JD_UTC, double* Ome_moon, double* Delta_Psi_1980,
                   double* epsilon_1980, double* epsilon_bar_1980); 

      Eigen::Matrix3f getNutationMat(double JD_UTC);

      Eigen::Matrix3f getPolarWMat(double JD_UTC);

      Eigen::MatrixXf getQmatrix(double deltaT);

      Eigen::Matrix3f getPrecessionMat(double JD_UTC);

      double getAlpha_G(double JD_UTC);

      Eigen::Matrix3f getS_primeMat(double JD_UTC);

      Eigen::Matrix3f getS_prime_dotMat(double JD_UTC);    
    
      Eigen::Matrix3f getECEF2ECI(double JD_UTC); // get matrix to convert from Cartesian inertial to Cartesian rotating  

      Eigen::Matrix3f getECEF2ECI_dot(double JD_UTC); // get derivative of matrix to convert from Cartesian inertial to Cartesian rotating  

      Eigen::Vector3f MoonPosition(double JD_TDB); // provides position of the Moon in inertial coordinates

      Eigen::Vector3f SunPosition(double JD_UT1, double JD_TDB); // provides position of the Sun in inertial coordinates
    
    private:

      void iCart2Kepl_(); // convert from Cartesian inertial state to Keplerian state

      Eigen::Vector3f Cart2spher_vec_(Eigen::Vector3f cartvec); // converts a vector from cartesian to spherical coordinates

      void Kepl2iCart_(); // convert from Keplerian state to Cartesian inertial state

      Eigen::Matrix3f Rot_axis_(double angle, int ax); // generate rotation matrix of angle around axis

      Eigen::Matrix3f orthoDCM_(Eigen::Matrix3f DCM); // orthogonalize matrix
 
      bool EXP_DRAGCANNON_ = false; // Flag for exponential cannonball density drag

      bool SRPCANNON_ = false; // Flag for cannonball srp

      bool LEG_GRAV_FIELD_ = false; // Flag for Legendre polynomial gravity field

      bool LUNISOLAR_ = false; // Flag for LuniSolar perturbations

      bool STM_PROP_ = false; // Flag for STM propagation
  };

}; 

#endif