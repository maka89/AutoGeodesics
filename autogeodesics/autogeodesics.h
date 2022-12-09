#pragma once

#include <array>

#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>

#include <Eigen/Eigen>
#include <tuple>
using namespace autodiff;
using namespace Eigen;

const double c_c{ 299792458.0 };
const double c_g{ 6.67430e-11 };


class AutoGeodesicsBase
{
protected:



    static VectorXvar matrix_to_vector(const Matrix4var& x);
    static Matrix4var vector_to_matrix(const VectorXvar& x);

    static VectorXvar metricfn2(const Vector4var& x, Matrix4var(*fn)(const Vector4var&));
    Matrix4var(*metfn_ptr)(const Vector4var& x);
   
    //NEWSTUFF
    std::tuple<std::array<Matrix4d, 4>, std::array<std::array<Matrix4d, 4>, 4> > calculate_christoffeld(const Vector4var& x);
    std::array<Matrix4d, 4> calculate_christoffel(const Vector4var& x);

    void get_christoffel(std::array<Matrix4d, 4>& chri_o, const VectorXd& vmetric, const MatrixXd& vjacobian);
    void get_christoffel(std::array<Matrix4d, 4>& chri_o, std::array<std::array<Matrix4d, 4>, 4>& dchri_o, const VectorXd& vmetric, const MatrixXd& vjacobian, const std::array<std::array<VectorXd, 4>, 4> &dvjac);
    Vector4d get_acc(const std::array<Matrix4d, 4>& chr, const Vector4d& velocity);
    
    Matrix4d vector_to_matrixd(const VectorXd& x);


    bool proper_time;

public:
    AutoGeodesicsBase() :metfn_ptr(NULL),proper_time(true) {};
    AutoGeodesicsBase(bool proper_time) :metfn_ptr(NULL), proper_time(proper_time) { };
    Vector4d calculate_acc(const Vector4d& x, const Vector4d& v);
    std::tuple<Vector4d, Matrix<double, 4, 8>> calculate_acc_jac(const Vector<double, 8>& inp);
    Vector4d setup_fourvelocity(const Vector4d &x, const Vector3d &v);

    void setMetFn(Matrix4var(*metfn_ptr)(const Vector4var& x)) { this->metfn_ptr = metfn_ptr; };

    
    
    

};


class AutoGeodesics : public AutoGeodesicsBase{
public:
    using AutoGeodesicsBase :: AutoGeodesicsBase;
    std::tuple<double, int> step_implicit_midpoint(Vector4d& x, Vector4d &v,  const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt, const double tol=1e-8);
    void step_rk4(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt);
    class Metrics;

    

private:
    void implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& now, const VectorXd& old, const double& dt);
    Vector<double, 8> rk_f(const Vector<double, 8>& d);
};

class AutoGeodesics::Metrics {
public:
    static inline Matrix4var schwarzschild_cartesian(const Vector4var& xi, double mass, Vector3d x0) {


        Vector4var x;
        x << xi[0], xi[1] - x0[0], xi[2] - x0[1], xi[3] - x0[2];

        var r2xy = x[1] * x[1] + x[2] * x[2];

        var r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
        var r = sqrt(r2);


        var A = 1.0 - 2.0 * c_g * mass / r / c_c / c_c;
        Matrix4var metric;


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
    static inline Matrix4var kerr_metric(const Vector4var& x, double mass, double angular_mom) {
        double rs = 2.0 * c_g * mass / pow(c_c, 2);
        double a = angular_mom / mass / c_c;
        double a2 = a * a;
        var r2 = x[1] * x[1];
        var sin2 = sin(x[2]) * sin(x[2]);

        var sigma = r2 + pow(a * sin(x[2]), 2);
        var delta = r2 - rs * x[1] + a2;

        Matrix4var metric = Matrix4var::Zero();

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
    static inline Matrix4var schwarzschild(const Vector4var& x, double mass) {

        double rs = 2.0 * c_g * mass / pow(c_c, 2);
        Matrix4var metric = Matrix4var::Zero();
        metric(0, 0) = -(1.0 - rs / x[1]);
        metric(1, 1) = 1.0 / (1.0 - rs / x[1]);
        metric(2, 2) = x[1] * x[1];
        metric(3, 3) = x[1] * x[1] * sin(x[2]) * sin(x[2]);

        return -1 * metric;

    };
};
