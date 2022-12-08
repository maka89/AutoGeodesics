//////
// Calculates the path of a test-particle at in orbit around earth
//////


#include <iostream>
#include <autogeodesics.h>
#include <chrono>
#include <string>
using namespace std::chrono;
using namespace std;


inline Matrix4var schwarzschild_cart(const Vector4var& x) {
    double mass = 5.972e24;
    Vector3d x0 = Vector3d::Zero();

    return AutoGeodesics::Metrics::schwarzschild_cartesian(x, mass, x0);
}

string format_method(bool rk) {
    if (rk)
        return string("Runge-Kutta 4");
    else
        return string("Implicit Midpoint Rule");
}

int main(int argc,char *argv[]) {
    bool rk = true;
    int steps = 1000;
    if (argc > 1) {
    
        string s;
        for (size_t i = 1; i < argc; i++) {
            s = string(argv[i]);
            if (s.find("-midpoint") < s.length()) {
                rk = false;
            }
            if (s.find("-niter") < s.length()) {
                steps = atol(argv[i + 1]);
            }
        }
    
    }
    Vector4d x, velocity, velo, xo;

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
    
    double dt = 2 * 3.141592 * sqrt(x[2] / 9.81997) / steps;
    auto start = high_resolution_clock::now();
    
    for (int i = 0; i < steps; i++) {
        velo = velocity;
        xo = x;


        if(rk)
            ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);
        else
            auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-3);


        //std::cout << x[0] / c_c << ", " << x[1] << ", " << x[2] << ", " << x[3] << std::endl;

    }
    auto duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    std::cerr << "Method: " << format_method(rk) << std::endl;
    std::cerr << "Steps: " << steps << std::endl;
    std::cerr << "Duration: " << duration.count() * 1e-6 << std::endl;

}
