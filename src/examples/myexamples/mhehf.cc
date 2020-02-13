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
    real_function_3d Vpsi=(V*psi);
    Vpsi.scale(-2.0).truncate();//2*Vpsi
    real_convolution_3d Gmu = BSHOperator3D(world,sqrt(-2*eps),0.001,1e-6);//G(mu) mu = sqrt(-2*E)
    real_function_3d tmp = apply(Gmu,Vpsi).truncate();// 
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps-0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}


int main(int argc, char** argv){
    // Initialization and finalizing the MADNESS environment
    //  MADNESS has its own parallel runtime that is fully compatible with MPI (Distributed memory)
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(1);//I think this is truncate f and df less than thresh
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);// This sets the domain

    vector_real_function_3d Vnuc=real_factory_3d(world).f(V).truncate_mode(0);//truncate V less than thresh
    vector_real_function_3d psi=real_factory_3d(world).f(guess);// inital psi
    psi.scale(1.0/psi.norm2());//vector_real_functions have the norm
    real_convolution_3d op = CoulombOperator(world,0.001,1e-6);

    double eps = -1.0;
    for (int iter =0; iter < 10; iter++){
        real_function_3d rho = square(psi).truncate();// what does truncate do here?
        real_function_3d potential = Vnuc + apply(op,rho).truncate();//what is apply
        iterate(world,potential,psi,eps);
    }

    double kinetic_energy = 0.0;
    for (int axis =0; axis < 3; axis++){
        real_derivative_3d D = free_space_derivative<double,3>(world,axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }

    real_funcimpl_3d rho = square(psi).truncate();
    double two_electron_energy = inner(apply(op,rho),rho);
    double nuclear_attraction_energy = 2.0*inner(Vnuc*psi,psi);
    double total_energy= kinetic_energy+nuclear_attraction_energy+two_electron_energy:
    
     // Manually tabluate the orbital along a line ... probably easier
    // to use the lineplot routine
    coord_3d r(0.0);
    psi.reconstruct();
    for (int i=0; i<201; i++) {
        r[2] = -L/2 + L*i/200.0;
        print(r[2], psi(r));
    }

   if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", ((nuclear_attraction_energy + two_electron_energy) / kinetic_energy));
    }

    world.gop.fence();
    finalize();
    return 0;

}