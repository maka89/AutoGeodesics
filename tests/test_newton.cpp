// C++ includes
#include <iostream>
#include <autogeodesics.h>
#include <vector>
#include <array>
#include <chrono>
#include <newton.hpp>

const double c_mass = 5.972e24;
const Vector3d c_x0 = { 0.0,0.0,0.0 };
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



std::vector<std::array<double, 6>> run_gr(int steps) {
    Vector4d xo, velo;

    AutoGeodesics ag = AutoGeodesics(false);

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

    Vector4d x = { 0.0,0.0,6371000.0,0.0 };
    double v = sqrt(9.81997 * x[2]);
    Vector4d velocity = ag.setup_fourvelocity(x, Vector3d({ v,0.0,0.0 })); //Turn velocity into 4-velocity.

    double t_end = 2 * 3.14159265 * sqrt(x[2] / 9.81997);
    double dt = t_end / steps;

    double err = -1.0;
    std::vector<std::array<double, 6>> data;

    for (int i = 0; i < steps; i++) {
        velo = velocity;
        xo = x;

        auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-9);
        //std::cout << err << ", " << niter << std::endl;
        //ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);

        std::array<double, 6> tmp;
        tmp = { (double)i + 1,(i + 1) * dt,x[0] / c_c,x[1] ,x[2],x[3] };
        data.push_back(tmp);

    }


    return data;
}

std::vector<std::array<double, 6>> run_newton(int steps) {
    Vector3d xo, velo;

    Newton ag = Newton();
    Vector3d x = { 0.0,6371000.0,0.0 };

    double v = sqrt(9.81997 * x[1]);
    ag.add_mass(c_mass, Vector3d({ 0,0,0 }));
   
    Vector3d velocity = { v,0.0,0.0 };

    double t_end = 2 * 3.14159265 * sqrt(x[1] / 9.81997);
    double dt = t_end / steps;


    double err = -1.0;
    std::vector<std::array<double, 6>> data;
    for (int i = 0; i < steps; i++) {
        velo = velocity;
        xo = x;

        auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-9);
        //std::cout << "Newton: "<< err << std::endl;

        //ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);
        std::array<double, 6> tmp;
        tmp = { (double)i + 1,(i + 1) * dt,(i + 1) * dt,x[0] ,x[1],x[2] };
        data.push_back(tmp);

    }
    

    return data;
}

int main() {
    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    std::vector<std::array<double, 6>> data1 = run_gr(20000);
    auto duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    std::cerr << "GR duration: " << duration.count() * 1e-6 << std::endl;

    start = high_resolution_clock::now();
    std::vector<std::array<double, 6>> data2 = run_newton(20000);
    duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    std::cerr << "Newton duration: " << duration.count() * 1e-6 << std::endl;


    double maxer = -1.0;
    double maxert = -1.0;
    for (size_t i = 0; i < data1.size(); i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < 3; j++)
            tmp += pow(data1[i][3 + j] - data2[i][3 + j], 2);
        tmp = sqrt(tmp);
        if (tmp > maxer)
            maxer = tmp;
        if (abs(data1[i][2] - data2[i][2]) > maxert)
            maxert = abs(data1[i][2] - data2[i][2]);

    }

    std::cout << data1[data1.size() - 1][3] << ", " << data1[data1.size() - 1][4] << ", " << data1[data1.size() - 1][5] << std::endl;
    std::cout << data2[data2.size() - 1][3] << ", " << data2[data2.size() - 1][4] << ", " << data2[data2.size() - 1][5] << std::endl;
    std::cout << "max error" << maxer << std::endl;
    std::cout << "max error time" << maxert << std::endl;
}
