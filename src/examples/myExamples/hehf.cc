#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

using namespace madness;


static const double L = 32.0;// box size
static const long k = 8;// box size
static const double thresh = 1e-6;// box size

//The intial guess exp(-2r) r=sqrt(x^2+y^2+z^2)
/* Q1 What is the importance of 1e-4
   Q2 What is the coord_3d
*/
static double guess(const coord_3d& r){
    const double x=r[0], y = r[1], z= r[2];
    return 6.0*exp(-2.0*sqrt(x*x+y*y+z*z+1e-4));
}

static double V(const coord_3d& r) {
    const double x=r[0], y = r[1], z= r[2];
    return -2.0/(sqrt(x*x+y*y+z*z+1e-8));
}


void iterate(World& world, real_function_3d V, real_function_3d& psi, double& eps){
    
}