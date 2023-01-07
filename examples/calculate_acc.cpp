// C++ includes
#include <iostream>
#include <autogeodesics.h>


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


int main() {
    Vector4d x,velocity;

    // Setup. proper_time=false. Use coordinate time.
    AutoGeodesics ag = AutoGeodesics(false);

    //Set metric function.
    ag.setMetFnComp(schwarzschild_cartesian_00, 0);
    //ag.setMetFnComp(schwarzschild_cartesian_01, 1);
    //ag.setMetFnComp(schwarzschild_cartesian_02, 2);
    //ag.setMetFnComp(schwarzschild_cartesian_03, 3);
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
    velocity = ag.setup_fourvelocity(x, Vector3d({ 0.0,0.0,0.0 })); 


    //Calculates acceleration. (0 component included).
    std::cout << ag.calculate_acc(x, velocity).cast<double>() << std::endl; 

}
