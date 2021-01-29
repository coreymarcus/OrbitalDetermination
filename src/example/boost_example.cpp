#include <iostream>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

const double gam = 0.15;

/* The rhs of x' = f(x) */
void lorenz( const state_type &x , state_type &dx ,  double t )
{
    dx[0] =  x[1];
    dx[1] = -x[0] - gam*x[1];
}


int main(int argc, char **argv)
{
    const double dt = 0.1;
    runge_kutta_dopri5<state_type> stepper;
    state_type x(2);
    x[0] = 1.0;
    x[1] = 0.0;


   double t = 0.0;
   cout << x[0] << endl;
   for ( size_t i(0); i <= 100; ++i){
       stepper.do_step(lorenz, x , t, dt );
       t += dt;
       cout << x[0] << endl;
   }


    return 0;
}