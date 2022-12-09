//////
// Calculates the path of a test-particle at in orbit around earth
//////


#include <iostream>
#include <autogeodesics.h>


inline Matrix4var schwarzschild_cart(const Vector4var& x) {
    double mass = 5.972e24;
    Vector3d x0 = Vector3d::Zero();

    return AutoGeodesics::Metrics::schwarzschild_cartesian(x, mass, x0);
}


int main() {
    Vector4d x, velocity,velo,xo;

    // Setup. proper_time=false. Use coordinate time.
    AutoGeodesics ag = AutoGeodesics(false);

    //Set metric function.
    ag.setMetFn(schwarzschild_cart);

    // Position y=Radius of Earth.
    x = { 0.0,0.0,6371000.0,0.0 };

    // Function to setup all 4 components of velocity. Mostly useful for proper_time=true.
    // In this case velocity could be set up with velocity << c,0,0,0;
    velocity = ag.setup_fourvelocity(x, Vector3d({ sqrt(9.81997 * x[2]),0.0,0.0 }));


    
    //Calculate 1 circular orbit at earth radius. Use 1000 steps.
    int steps = 1000;
    double dt= 2 * 3.141592 * sqrt(x[2] / 9.81997) / steps;

    for (int i = 0; i < steps; i++) {
        velo = velocity;
        xo = x;

        
        // Make one step using RungeKutta4
        //ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);


        //Alternatively, use implicit midpoint rule
        auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-6);


        std::cout << x[0]/c_c<<", "<<x[1]<<", "<<x[2]<<", "<<x[3] << std::endl;

    }


}
