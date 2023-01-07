//////
// Calculates the path of a test-particle at in orbit around earth
//////


#include <iostream>
#include <autogeodesics.h>
#include <chrono>
#include <string>
using namespace std::chrono;
using namespace std;


inline dual2nd schwarzschild_cartesian_00(Vector4dual2nd& xi) { double mass = 5.972e24;  Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return A; }
inline dual2nd schwarzschild_cartesian_01(Vector4dual2nd& xi) { return 0; }
inline dual2nd schwarzschild_cartesian_02(Vector4dual2nd& xi) { return 0; }
inline dual2nd schwarzschild_cartesian_03(Vector4dual2nd& xi) { return 0; }
inline dual2nd schwarzschild_cartesian_11(Vector4dual2nd& xi) { double mass = 5.972e24; Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return -1.0 * (x[1] * x[1] / A / r2 + (1.0 / r2) * x[1] * x[1] * x[3] * x[3] + x[2] * x[2] / r2xy); }
inline dual2nd schwarzschild_cartesian_12(Vector4dual2nd& xi) { double mass = 5.972e24; Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return -1.0 * 0.5 * (2.0 * (x[1] * x[2] / r2) / A + (1.0 / r2 / r2xy) * 2.0 * x[1] * x[2] * x[3] * x[3] + -2.0 * x[1] * x[2] / r2xy); }
inline dual2nd schwarzschild_cartesian_13(Vector4dual2nd& xi) { double mass = 5.972e24; Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return -1.0 * 0.5 * (2.0 * (x[1] * x[3] / r2) / A + -2.0 * x[1] * x[3] / r2); }
inline dual2nd schwarzschild_cartesian_22(Vector4dual2nd& xi) { double mass = 5.972e24; Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return -1.0 * (x[2] * x[2] / A / r2 + (1.0 / r2) * x[2] * x[2] * x[3] * x[3] + x[1] * x[1] / r2xy); }
inline dual2nd schwarzschild_cartesian_23(Vector4dual2nd& xi) { double mass = 5.972e24; Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return -1.0 * 0.5 * (2.0 * (x[2] * x[3] / r2) / A + -2.0 * x[2] * x[3] / r2); }
inline dual2nd schwarzschild_cartesian_33(Vector4dual2nd& xi) { double mass = 5.972e24; Vector4dual2nd x; x << xi[0], xi[1], xi[2], xi[3]; dual2nd r2xy = x[1] * x[1] + x[2] * x[2]; dual2nd r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3]; dual2nd r = sqrt(r2); dual2nd A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c; return -1.0 * (x[3] * x[3] / A / r2 + (r2xy / r2)); }



inline Matrix4dual2nd schwarzschild_cart(const Vector4dual2nd& x) {
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
    int steps = 5000;
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
    ag.setMetFnComp(schwarzschild_cartesian_00, 0);
    //ag.setMetFnComp(schwarzschild_cartesian_01, 1);
    //ag.setMetFnComp(schwarzschild_cartesian_02, 2);
    //ag.setMetFnComp(schwarzschild_cartesian_03,  3);
    ag.setMetFnComp(schwarzschild_cartesian_11, 4);
    ag.setMetFnComp(schwarzschild_cartesian_12, 5);
    ag.setMetFnComp(schwarzschild_cartesian_13, 6);
    ag.setMetFnComp(schwarzschild_cartesian_22, 7);
    ag.setMetFnComp(schwarzschild_cartesian_23, 8);
    ag.setMetFnComp(schwarzschild_cartesian_33, 9);

    // Position y=Radius of Earth.
    x = { 0.0,0.0,6371000.0,0.0 };

    // Function to setup all 4 components of velocity. Mostly useful for proper_time=true.
    // In this case velocity could be set up with velocity << c,0,0,0;
    std::cerr << "Method 1" << std::endl;
    velocity = ag.setup_fourvelocity(x, Vector3d({ sqrt(9.81997 * x[2]),0.0,0.0 }));
    std::cerr << "Method 2" << std::endl;


    //Calculate 1 circular orbit at earth radius. Use 1000 steps.
    
    double dt = 2 * 3.141592 * sqrt(x[2] / 9.81997) / steps;
    auto start = high_resolution_clock::now();
    
    for (int i = 0; i < steps; i++) {
        velo = velocity;
        xo = x;


        if(rk)
            ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);
        else
            auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-6);


        //std::cout << x[0] / c_c << ", " << x[1] << ", " << x[2] << ", " << x[3] << std::endl;

    }
    auto duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    std::cerr << "Method: " << format_method(rk) << std::endl;
    std::cerr << "Steps: " << steps << std::endl;
    std::cerr << "Duration: " << duration.count() * 1e-6 << std::endl;

}
