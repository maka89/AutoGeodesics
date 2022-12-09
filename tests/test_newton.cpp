// C++ includes
#include <iostream>
#include <autogeodesics.h>
#include <vector>
#include <array>
#include <chrono>
#include <newton.hpp>

const double c_mass = 5.972e24;
const Vector3d c_x0 = { 0.0,0.0,0.0 };

inline var potential_fn(const Vector4var& x) {
    Vector3d x0({ 0.0,0.0,0.0 });
    var r2 = pow(x[1] - x0[0], 2) + pow(x[2] - x0[1], 2) + pow(x[3] - x0[2], 2);
    var r= sqrt(r2);
    return -c_mass / r;
}
inline Matrix4var schwarzschild_cart(const Vector4var &x) {
    double mass = c_mass;
    Vector3d x0 = c_x0;

    return AutoGeodesics::Metrics::schwarzschild_cartesian(x, mass, x0);
}


std::vector<std::array<double,6>> run_gr(int steps) {
    Vector4d xo, velo;
    
    AutoGeodesics ag = AutoGeodesics(false);
    
    ag.setMetFn(schwarzschild_cart);
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

        auto [err,niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo,velo), dt, 1e-3);
        //std::cout << err << ", " << niter << std::endl;
        //ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);

        std::array<double, 6> tmp;
        tmp = { (double)i+1,(i+1) * dt,x[0]/c_c,x[1] ,x[2],x[3] };
        data.push_back(tmp);

    }


    return data;
}

std::vector<std::array<double, 6>> run_newton(int steps) {
    Vector4d xo, velo;

    Newton ag = Newton();
    Vector4d x = { 0.0,0.0,6371000.0,0.0 };
    double v = sqrt(9.81997 * x[2]);
    ag.set_potential(&potential_fn);

    Vector4d velocity = { 0.0,v,0.0,0.0 };

    double t_end = 2 * 3.14159265 * sqrt(x[2] / 9.81997);
    double dt = t_end / steps;

    double err = -1.0;
    std::vector<std::array<double, 6>> data;

    for (int i = 0; i < steps; i++) {
        velo = velocity;
        xo = x;

        //auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-9);
        ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);
        std::array<double, 6> tmp;
        tmp = { (double)i + 1,(i + 1) * dt,(i + 1) * dt,x[1] ,x[2],x[3] };
        data.push_back(tmp);

    }


    return data;
}

int main() {
    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    std::vector<std::array<double, 6>> data1 = run_gr(200000);
    auto duration = duration_cast<microseconds>(high_resolution_clock::now() - start);
    std::cerr << "GR duration: " << duration.count() * 1e-6 << std::endl;

    start = high_resolution_clock::now();
    std::vector<std::array<double, 6>> data2 = run_newton(200000);
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
    std::cout << "max error" << maxer<<std::endl;
    std::cout << "max error time" << maxert << std::endl;
}
