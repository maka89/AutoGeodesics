#pragma once

#include <array>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <Eigen/Eigen>
#include <tuple>
using namespace autodiff;
using namespace Eigen;

const double c_c{ 299792458.0 };
const double c_g{ 6.67430e-11 };


class AutoGeodesicsBase
{
protected:
    Vector4d get_acc(std::array<Matrix4d, 4> chr, const Vector4d& velocity);
    Vector4real get_acc(std::array<Matrix4d, 4> chr, const Vector4real& velocity);

    std::array<Matrix4d, 4> get_christoffel(VectorXd vmetric, MatrixXd vjacobian);

    static VectorXreal matrix_to_vector(const Matrix4real& x);
    static Matrix4d vector_to_matrix(const VectorXd& x);

    static VectorXreal metricfn2(const Vector4real& x, Matrix4real(*fn)(const Vector4real&));
    std::array<Matrix4d, 4> calculate_christoffel(const Vector4d& x);
    Matrix4real(*metfn_ptr)(const Vector4real& x);
public:
    AutoGeodesicsBase() :metfn_ptr(NULL) {};

    Vector4d setup_fourvelocity(const Vector4d &x, const Vector3d &v, bool input_is_four_velocity = false);
    Vector4d calculate_acc(const Vector4d &x, const Vector4d &velocity );

    void setMetFn(Matrix4real(*metfn_ptr)(const Vector4real& x)) { this->metfn_ptr = metfn_ptr; };

    Vector4d calculate_acc_newton(const Vector3d &x, const double &mass, const Vector3d &x0); 
    class Metrics;

};


class AutoGeodesics : public AutoGeodesicsBase{
public:
    void step_euler_midpoint(Vector4d& x, Vector4d& v, Vector4d& acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old, const double& dt);
    void step_euler(Vector4d& x, Vector4d& v, Vector4d& acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old,  const double& dt);
    void step_euler_cromer(Vector4d& x, Vector4d& v, Vector4d& acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old,  const double& dt);
    double step_velocity_verlet(Vector4d& x, Vector4d &v, Vector4d &acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old, const double& dt, const double tol=1e-8);

    void calculate_acc_velocity_jacobian(Vector4d& acc, Matrix4d& J, const Vector4d& velocity, const Vector4d& x); 
    
private:
    
    int velocity_verlet_resids(Vector4d& resids, Matrix4d& jac, const Vector4d& velocity, const std::array<Matrix4d, 4>& chr, const Vector4d& v_old, const Vector4d& acc_old, const Vector4d& x, const double& dt);
};

class AutoGeodesicsBase::Metrics {
public:
    static inline Matrix4real schwarzschild_cartesian(const Vector4real& xi, double mass, Vector3d x0) {


        Vector4real x;
        x << xi[0], xi[1] - x0[0], xi[2] - x0[1], xi[3] - x0[2];

        Real r2xy = x[1] * x[1] + x[2] * x[2];

        Real r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
        Real r = sqrt(r2);


        Real A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c;
        Matrix4real metric;


        for (size_t i = 0; i < 3; i++) {
            metric(0, i + 1) = 0.0;
            metric(i + 1, 0) = 0.0;
        }

        metric(0, 0) = -A;

        metric(1, 1) = x[1] * x[1] / A / r2;
        metric(2, 2) = x[2] * x[2] / A / r2;
        metric(3, 3) = x[3] * x[3] / A / r2;
        metric(1, 2) = 2.0 * (x[1] * x[2] / r2) / A;
        metric(1, 3) = 2.0 * (x[1] * x[3] / r2) / A;
        metric(2, 3) = 2.0 * (x[2] * x[3] / r2) / A;

        metric(1, 1) += (1.0 / r2) * x[1] * x[1] * x[3] * x[3];
        metric(2, 2) += (1.0 / r2) * x[2] * x[2] * x[3] * x[3];
        metric(3, 3) += (r2xy / r2);
        metric(1, 2) += (1.0 / r2 / r2xy) * 2.0 * x[1] * x[2] * x[3] * x[3];
        metric(1, 3) += -2.0 * x[1] * x[3] / r2;
        metric(2, 3) += -2.0 * x[2] * x[3] / r2;

        metric(1, 1) += x[2] * x[2] / r2xy;
        metric(2, 2) += x[1] * x[1] / r2xy;
        metric(1, 2) += -2.0 * x[1] * x[2] / r2xy;


        metric(1, 2) *= 0.5;
        metric(1, 3) *= 0.5;
        metric(2, 3) *= 0.5;

        metric(2, 1) = 1.0 * metric(1, 2);
        metric(3, 1) = 1.0 * metric(1, 3);
        metric(3, 2) = 1.0 * metric(2, 3);


        return -1.0 * metric;

    };
    static inline Matrix4real kerr_metric(const Vector4real& x, double mass, double angular_mom) {
        double rs = 2.0 * c_g * mass / pow(c_c, 2);
        double a = angular_mom / mass / c_c;
        double a2 = a * a;
        Real r2 = x[1] * x[1];
        Real sin2 = sin(x[2]) * sin(x[2]);

        Real sigma = r2 + pow(a * sin(x[2]), 2);
        Real delta = r2 - rs * x[1] + a2;

        Matrix4real metric = Matrix4real::Zero();

        metric(0, 0) = -(1.0 - (rs * x[1]) / sigma);
        metric(1, 1) = sigma / delta;
        metric(2, 2) = sigma;
        metric(3, 3) = (r2 + a2 + rs * x[1] * a2 * sin2 / sigma) * sin2;


        metric(0, 3) = -2.0 * rs * x[1] * a * sin2 / sigma;

        //DIVIDE BY 2? Line element vs metric... Off diagonal elements.. Same for schwarzchild cart?
        metric(0, 3) *= 0.5;

        metric(3, 0) = metric(0, 3);

        return -1.0 * metric;
    };
    static inline Matrix4real schwarzschild(const Vector4real& x, double mass) {

        double rs = 2.0 * c_g * mass / pow(c_c, 2);
        Matrix4real metric = Matrix4real::Zero();
        metric(0, 0) = -(1.0 - rs / x[1]);
        metric(1, 1) = 1.0 / (1.0 - rs / x[1]);
        metric(2, 2) = x[1] * x[1];
        metric(3, 3) = x[1] * x[1] * sin(x[2]) * sin(x[2]);

        return -1 * metric;

    };
};
