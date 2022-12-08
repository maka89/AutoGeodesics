#pragma once
#include "autogeodesics.h"
#include <iostream>
inline VectorXvar AutoGeodesicsBase::metricfn2(const Vector4var& x,Matrix4var(*fn)(const Vector4var &)) {
    return matrix_to_vector(fn(x));
}


inline VectorXvar AutoGeodesicsBase::matrix_to_vector(const Matrix4var& x) {
    VectorXvar v(16);

    size_t  k = 0;
    for (size_t  i = 0; i < 4; i++) {
        for (size_t  j = 0; j < 4; j++) {
            v[k] = x(i, j);
            k++;
        }
    }
    return v;
}



Matrix4var AutoGeodesicsBase::vector_to_matrix(const VectorXvar& x) {
    Matrix4var m;

    size_t  k = 0;
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            m(i, j) = x[k];
            k++;
        }
    }
    

    return m;
}

std::array<Matrix4var, 4> AutoGeodesicsBase::get_christoffel(const VectorXvar &vmetric, const MatrixXvar &vjacobian) {
    Matrix4var metric;
    metric = vector_to_matrix(vmetric);



    std::array<Matrix4var, 4> dmetric;
    for (size_t i = 0; i < 4; i++)
        dmetric[i] = vector_to_matrix(vjacobian.col(i));
    


    std::array<Matrix4var, 4> christoffel;
    Vector4var xx;
    Matrix4var metinv = metric.inverse();

    
    for (size_t i = 0; i < 4; i++) {
        for (size_t k = 0; k < 4; k++) {
            for (size_t l = 0; l < 4; l++) {

                var tot = 0.0;
                for (size_t m = 0; m < 4; m++) {
                    tot += metinv(i, m) * (dmetric[l](m, k) + dmetric[k](m, l) - dmetric[m](k, l));
                }
                tot *= 0.5;
                christoffel[i](k, l) = tot;
            }
        }
    }

    
    return christoffel;
}


var AutoGeodesicsBase::dotdot(const Matrix4var &m, const Vector4var &v,bool symmetric) {

    var tmp = 0.0;
    if (symmetric) {
        for (size_t i = 0; i < 4; i++)
            tmp += m(i, i) * v[i] * v[i];
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = i + 1; j < 4; j++)
                tmp += 2.0 * m(i, j) * v[i] * v[j];
        }
        return tmp;
    }
    else
        return v.transpose() * m * v;
   
}

Vector4var AutoGeodesicsBase::get_acc(const std::array<Matrix4var, 4> &chr,const Vector4var &velocity) {
    Vector4var acc;
    acc = Vector4var::Zero();
    var tmp = dotdot(chr[0], velocity);
    acc[0] = -tmp;
    for (size_t i = 1; i < 4; i++)
        acc[i] = -dotdot(chr[i], velocity);

    if (!this->proper_time) {
        for (size_t i = 0; i < 4; i++)
            acc[i] += tmp * (velocity[i]/c_c);
    }
    return acc;

}


Vector4d AutoGeodesicsBase::setup_fourvelocity(const Vector4d& x, const Vector3d& velocity) {
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
    /*
    if (input_is_four_velocity) {
        double a = metric(0, 0);
        double b = c0sum;
        double c = vsum - 1.0;

        double lambda = 0.5 * (-b + sqrt(b * b - 4.0 * a * c)) / a;
        
        if (this->proper_time) {
            vv << lambda * c_c, velocity[0], velocity[1], velocity[2];
        }
        else {
            vv << c_c, velocity[0] * lambda, velocity[1] * lambda, velocity[2] * lambda;
        }
    }*/

   
    vv << c_c, velocity[0], velocity[1], velocity[2];
    if (this->proper_time) {
        double f = metric(0, 0);
        f += c0sum;
        f += vsum;
        double lambda = sqrt(1.0 / f);
        vv = lambda * vv;
    }



    
    
    return vv;
}
Vector4d AutoGeodesicsBase::calculate_acc(const Vector4d& x, const Vector4d& velocity) {

    Vector4var xx;
    xx << x[0], x[1], x[2], x[3];
    return calculate_acc(xx, velocity.cast<var>()).cast<double>();

}

Vector4var AutoGeodesicsBase::calculate_acc(const Vector4var& x, const Vector4var& velocity) {


    std::array<Matrix4var, 4> christoffel = calculate_christoffel(x);
    Vector4var acc4 = get_acc(christoffel, velocity);
    return acc4;

}

std::array<Matrix4var, 4> AutoGeodesicsBase::calculate_christoffel(const Vector4var& x) {
    VectorXvar vmetric = metricfn2(x,this->metfn_ptr);
    Matrix<var, 16, 4> J = Matrix<var, 16, 4>::Zero();
    for (size_t i = 0; i < 16; i++) {
        auto tmp = derivatives(vmetric(i), wrt(x[0],x[1],x[2],x[3]));
        for(size_t j=0;j<4;j++)
            J(i, j) = tmp[j];
    }

    std::array<Matrix4var, 4> christoffel;
    christoffel = get_christoffel(vmetric, J);
    return christoffel;


}
