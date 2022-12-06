#pragma once
#include "autogeodesics.h"

inline VectorXreal AutoGeodesicsBase::metricfn2(const Vector4real& x,Matrix4real(*fn)(const Vector4real &)) {
    return matrix_to_vector(fn(x));
}
Vector4d AutoGeodesicsBase::calculate_acc_newton(const Vector3d& x, const double& mass, const Vector3d& x0) {
    Vector3d dirvec = x - x0;
    double r2 = dirvec.transpose() * dirvec;
    dirvec /= sqrt(r2);

    Vector4d retval;
    retval[0] = 0.0;
    for (size_t i = 0; i < 3; i++) 
        retval[i+1] = -dirvec[i] * mass * c_g / r2;
    return retval;
}

inline VectorXreal AutoGeodesicsBase::matrix_to_vector(const Matrix4real& x) {
    VectorXreal v(16);

    size_t  k = 0;
    for (size_t  i = 0; i < 4; i++) {
        for (size_t  j = 0; j < 4; j++) {
            v[k] = x(i, j);
            k++;
        }
    }
    return v;
}

Matrix4d AutoGeodesicsBase::vector_to_matrix(const Eigen::VectorXd& x) {
    Eigen::Matrix4d m;

    size_t  k = 0;
    for (size_t  i = 0; i < 4; i++) {
        for (size_t  j = 0; j < 4; j++) {
            m(i, j) = x[k];
            k++;
        }
    }
    return m;
}

std::array<Matrix4d, 4> AutoGeodesicsBase::get_christoffel(VectorXd vmetric, MatrixXd vjacobian) {
    Matrix4d metric;
    metric = vector_to_matrix(vmetric);



    std::array<Matrix4d, 4> dmetric;
    for (size_t  i = 0; i < 4; i++)
        dmetric[i] = vector_to_matrix(vjacobian.col(i));

    std::array<Matrix4d, 4> christoffel;
    Vector4d xx;
    Matrix4d metinv = metric.inverse();


    for (size_t  i = 0; i < 4; i++) {
        for (size_t  k = 0; k < 4; k++) {
            for (size_t  l = 0; l < 4; l++) {

                double tot = 0.0;
                for (size_t  m = 0; m < 4; m++) {
                    tot += metinv(i, m) * (dmetric[l](m, k) + dmetric[k](m, l) - dmetric[m](k, l));
                }
                tot *= 0.5;
                christoffel[i](k, l) = tot;
            }
        }
    }
    return christoffel;
}

Vector4d AutoGeodesicsBase::get_acc(std::array<Matrix4d, 4> chr,const Vector4d &velocity ) {
    Vector4d acc;
    acc = Vector4d::Zero();
    for (size_t i = 0; i < 4; i++)
        acc[i] = -velocity.transpose() * chr[i] * velocity;
    return acc;

}

Vector4real AutoGeodesicsBase::get_acc(std::array<Matrix4d, 4> chr,const Vector4real &velocity) {
    Vector4real acc;
    acc = Vector4real::Zero();
    for (size_t i = 0; i < 4; i++) 
        acc[i] = -velocity.transpose()*chr[i] * velocity;
    return acc;

}


Vector4d AutoGeodesicsBase::setup_fourvelocity(const Vector4d& x, const Vector3d& velocity, bool input_is_four_velocity) {
    Vector4d vv;
    Matrix4d metric = this->metfn_ptr(x).cast<double>();
    Matrix3d metric3 = metric(Eigen::seq(1, last), Eigen::seq(1, last));
    Vector3d col0 = metric.col(0)(Eigen::seq(1, last));
    Vector3d vc = velocity / c_c;

    double vsum = 0.0;
    double c0sum = 0.0;
    for (size_t i = 0; i < 3; i++) {
        c0sum += 2.0 * col0[i];
        for (size_t j = 0; j < 3; j++) {
            vsum += metric3(i, j) * vc(i) * vc(j);
        }
    }

    if (input_is_four_velocity) {
        double a = metric(0, 0);
        double b = c0sum;
        double c = vsum - 1.0;

        double lambda = 0.5 * (-b + sqrt(b * b - 4.0 * a * c)) / a;
        vv << lambda * c_c, velocity[0], velocity[1], velocity[2];
    }
    else {
        double f = metric(0, 0);
        f += c0sum;
        f += vsum;
        double lambda = sqrt(1.0 / f);
        vv <<  c_c,  velocity[0], velocity[1], velocity[2];
        vv = lambda * vv;
    }
    return vv;
}
Vector4d AutoGeodesicsBase::calculate_acc(const Vector4d &x, const Vector4d &velocity) {

    std::array<Matrix4d, 4> christoffel = calculate_christoffel(x);
    Vector4d acc4 = get_acc(christoffel, velocity);
   
    return acc4;

}


std::array<Matrix4d, 4> AutoGeodesicsBase::calculate_christoffel(const Vector4d& x) {
    Vector4real xx;

    xx << x[0], x[1], x[2], x[3];

    VectorXreal vmetric;
    MatrixXd J;
    Vector4d ders;

    J = jacobian(&metricfn2, wrt(xx), at(xx,this->metfn_ptr), vmetric);
    Matrix4d metric = vector_to_matrix(vmetric.cast<double>());

    Vector4d acc4;

    std::array<Matrix4d, 4> christoffel;
    christoffel = get_christoffel(vmetric.cast<double>(), J.cast<double>());
    return christoffel;

}
