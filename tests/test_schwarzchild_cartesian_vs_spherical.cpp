// C++ includes
#include <iostream>
#include <autogeodesics.h>
#include <vector>
#include <array>

inline Matrix4real schwarzchild_cart(const Vector4real &x) {
    double mass = 5.972e24;
    Vector3d x0 = { 0.0,0.0,0.0 };

    return AutoGeodesics::Metrics::schwarzchild_cartesian(x, mass, x0);
}

 inline Matrix4real schwarzchild(const Vector4real& x) {
    double mass = 5.972e24;

    return  AutoGeodesics::Metrics::schwarzchild(x, mass);

}
std::vector<std::array<double,5>> run( bool cart) {
    Vector4d x, vel, acc;
    Vector4d xo, velo, acco;
    Vector3d  vel3;

    

    Matrix4real(*metricfn)(const Vector4real & x);
    double escape_vel = 0.99*c_c;
    double r0 = 0.02;// 6371000.0;
    double m_pi = 3.1415926535897932384626433;
    if (cart) {
        metricfn = &schwarzchild_cart;
        x << 0.0, r0 / sqrt(2), r0 / sqrt(2), 0.0;
        vel3 << escape_vel/sqrt(2), escape_vel/sqrt(2), 0.0;
    }
    else {
        metricfn = &schwarzchild;
        double theta = 0.5 * m_pi;
        double phi = 0.25 * m_pi;
        x << 0.0, r0, theta, phi;
        vel3 << escape_vel, 0.0,0.0;
    }

    AutoGeodesics ag = AutoGeodesics();
    ag.setMetFn(metricfn);
    vel = ag.setup_fourvelocity(x, vel3, true);
    acc = ag.calculate_acc(x, vel);

    int steps = 20000;
    double t_end = 20.0/c_c;
    double dt = t_end / steps;

    double err = -1.0;
    std::vector<std::array<double, 5>> data;

    for (int i = 0; i < steps; i++) {
        acco = acc;
        velo = vel;
        xo = x;

        double tol = 1e-6;
        err = ag.step_velocity_verlet(x, vel, acc, std::make_tuple(xo, velo, acco), dt, tol);
      

        std::array<double, 5> tmp;
        if (cart)
            tmp = { (double)i,i * dt,x[1] ,x[2],x[3] };
        else
            tmp = { (double)i,i * dt,x[1] * sin(x[2]) * cos(x[3]) ,x[1] * sin(x[2]) * sin(x[3]),x[1] * cos(x[2]) };
        data.push_back(tmp);

    }


    return data;
}

int main() {
    std::vector<std::array<double, 5>> data1 = run(false);
    std::vector<std::array<double, 5>> data2 = run(true);

    double maxer = -1.0;
    for (size_t i = 0; i < data1.size(); i++) {
        double tmp = 0.0;
        for (size_t j = 0; j < 3; j++)
            tmp += pow((data1[i][2 + j] - data2[i][2 + j]), 2);
        tmp = sqrt(tmp);
        if (tmp > maxer)
            maxer = tmp;

    }
    std::cout << "max error" << maxer<<std::endl;
}