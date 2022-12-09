#pragma once
#include "autogeodesics.h"
#include <iostream>
inline VectorXvar AutoGeodesicsBase::metricfn2(const Vector4var& x,Matrix4var(*fn)(const Vector4var &)) {
    return matrix_to_vector(fn(x));
}

/*
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



Matrix4d AutoGeodesicsBase::vector_to_matrixd(const VectorXd& x) {
    Matrix4d m;

    size_t  k = 0;
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            m(i, j) = x[k];
            k++;
        }
    }


    return m;
}*/


inline VectorXvar AutoGeodesicsBase::matrix_to_vector(const Matrix4var& x) {
    VectorXvar v(10);

    size_t  k = 0;
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = i; j < 4; j++) {
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
        for (size_t j = i; j < 4; j++) {
            m(i, j) = x[k];
            m(j, i) = x[k];
            k++;
        }
    }


    return m;
}



Matrix4d AutoGeodesicsBase::vector_to_matrixd(const VectorXd& x) {
    Matrix4d m;

    size_t  k = 0;
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = i; j < 4; j++) {
            m(i, j) = x[k];
            m(j, i) = x[k];
            k++;
        }
    }


    return m;
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

void AutoGeodesicsBase::get_christoffel(std::array<Matrix4d, 4> &chri_o, const VectorXd& vmetric, const MatrixXd& vjacobian) {


    Matrix4d metric;
    metric = vector_to_matrixd(vmetric);



    std::array<Matrix4d, 4> dmetric;
    for (size_t i = 0; i < 4; i++)
        dmetric[i] = vector_to_matrixd(vjacobian.col(i));


    std::array<Matrix4d, 4> christoffel;
    Vector4d xx;
    Matrix4d metinv = metric.inverse();


    for (size_t i = 0; i < 4; i++) {
        for (size_t k = 0; k < 4; k++) {
            for (size_t l = 0; l < 4; l++) {

                double tot = 0.0;
                for (size_t m = 0; m < 4; m++) {
                    tot += metinv(i, m) * (dmetric[l](m, k) + dmetric[k](m, l) - dmetric[m](k, l));
                }
                tot *= 0.5;
                christoffel[i](k, l) = tot;
            }
        }
    }
    chri_o = christoffel;
    return;
}

void AutoGeodesicsBase::get_christoffel(std::array<Matrix4d, 4> &chri_o,std::array<std::array<Matrix4d, 4>, 4> &dchri_o,const VectorXd& vmetric, const MatrixXd& vjacobian, const std::array<std::array<VectorXd,4>,4> &dvjac) {
    Matrix4d metric;
    metric = vector_to_matrixd(vmetric);



    std::array<Matrix4d, 4> dmetric;
    for (size_t i = 0; i < 4; i++)
        dmetric[i] = vector_to_matrixd(vjacobian.col(i));

   
    std::array<Matrix4d, 4> christoffel;
    Vector4d xx;
    Matrix4d metinv = metric.inverse();


    for (size_t i = 0; i < 4; i++) {
        for (size_t k = 0; k < 4; k++) {
            for (size_t l = 0; l < 4; l++) {

                double tot = 0.0;
                for (size_t m = 0; m < 4; m++) {
                    tot += metinv(i, m) * (dmetric[l](m, k) + dmetric[k](m, l) - dmetric[m](k, l));
                }
                tot *= 0.5;
                christoffel[i](k, l) = tot;
            }
        }
    }
    chri_o = christoffel;

    // DERIVATIVE
        

    std::array<std::array<Matrix4d, 4>, 4> ddmetric;
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++)
            ddmetric[i][j] = vector_to_matrixd(dvjac[i][j]);
    }

    std::array<Matrix4d, 4> dmetinv;
    for (size_t i = 0; i < 4; i++)
        dmetinv[i] = -metinv * dmetric[i] * metinv;

    std::array<std::array<Matrix4d, 4>, 4> dchris;
    for (size_t u = 0; u < 4; u++) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t k = 0; k < 4; k++) {
                for (size_t l = 0; l < 4; l++) {

                    double tot = 0.0;
                    for (size_t m = 0; m < 4; m++) {
                        tot += dmetinv[u](i, m) * (dmetric[l](m, k) + dmetric[k](m, l) - dmetric[m](k, l));
                    }
                    tot *= 0.5;
                    dchris[u][i](k, l) = tot;

                    tot = 0.0;
                    for (size_t m = 0; m < 4; m++) {
                        tot += metinv(i, m) * (ddmetric[u][l](m, k) + ddmetric[u][k](m, l) - ddmetric[u][m](k, l));
                    }
                    tot *= 0.5;
                    dchris[u][i](k, l) += tot;

                }
            }
        }
    }

    dchri_o = dchris;
    
    return;
}

std::tuple<std::array<Matrix4d, 4>, std::array<std::array<Matrix4d, 4>, 4> > AutoGeodesicsBase::calculate_christoffeld(const Vector4var& xx) {
    Vector4var x = xx;
    VectorXvar vmetric = metricfn2(x, this->metfn_ptr);
    Matrix<double, 10, 4> J = Matrix<double, 10, 4>::Zero();
    std::array<std::array<VectorXd, 4>, 4> H;

    for (size_t i = 0; i < 4; i++)for (size_t j = 0; j < 4; j++)H[i][j] = Vector<double, 10>::Zero();

    for (size_t i = 0; i < 10; i++) {
        auto tmp = derivativesx(vmetric(i), wrt(x[0], x[1], x[2], x[3]));
        for (size_t j = 0; j < 4; j++)
            J(i, j) = (double)tmp[j];
        
        for (size_t k = 0; k < 4; k++) {
            auto h = derivatives(tmp[k], wrt(x[0], x[1], x[2], x[3]));
            for (size_t j = 0; j < 4; j++)
                H[k][j](i) = (double)h[j];
        }

    }

    std::array<std::array<Matrix4d, 4>, 4> dchris;
    std::array<Matrix4d, 4> christoffel;
    get_christoffel(christoffel,dchris,vmetric.cast<double>(), J,H);
    return std::make_tuple(christoffel,dchris);


}
std::array<Matrix4d, 4> AutoGeodesicsBase::calculate_christoffel(const Vector4var& x) {
    VectorXvar vmetric = metricfn2(x, this->metfn_ptr);
    Matrix<double, 10, 4> J = Matrix<double, 10, 4>::Zero();
    

    for (size_t i = 0; i < 10; i++) {
        auto tmp = derivativesx(vmetric(i), wrt(x[0], x[1], x[2], x[3]));
        for (size_t j = 0; j < 4; j++)
            J(i, j) = (double)tmp[j];
    }

    std::array<Matrix4d, 4> christoffel;
    get_christoffel(christoffel,  vmetric.cast<double>(), J);
    return christoffel;
}


Vector4d AutoGeodesicsBase::get_acc(const std::array<Matrix4d, 4>& chr, const Vector4d& velocity) {
    Vector4d acc;
    acc = Vector4d::Zero();
    double tmp = velocity.transpose() * chr[0] * velocity;
    acc[0] = -tmp;
    for (size_t i = 1; i < 4; i++)
        acc[i] = -velocity.transpose()*chr[i]*velocity;

    if (!this->proper_time) {
        for (size_t i = 0; i < 4; i++)
            acc[i] += tmp * (velocity[i] / c_c);
    }
    return acc;

}
Vector4d AutoGeodesicsBase::calculate_acc(const Vector4d& x, const Vector4d& v) {
    std::array<Matrix4d, 4> chr = calculate_christoffel(x);
    return get_acc(chr, v);
}
std::tuple<Vector4d,Matrix<double,4,8>> AutoGeodesicsBase::calculate_acc_jac( const Vector<double,8> &inp) {
    Vector4d x, velocity;
    x << inp[0], inp[1], inp[2], inp[3];
    velocity << inp[4], inp[5], inp[6], inp[7];

    auto [christoffel , dchris] = calculate_christoffeld(x);
    Vector4d acc4 = get_acc(christoffel, velocity);
    Matrix<double, 4, 8> jac;

    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            jac(i, j) = -velocity.transpose() * dchris[j][i] * velocity;
        }
    }

    if (!proper_time) {
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                double tmp = velocity.transpose() * dchris[j][0] * velocity;
                jac(i, j) += tmp* (velocity[i] / c_c);
            }
        }
    }

    //Velocity
    for (size_t i = 0; i < 4; i++) {
        auto tmp = -2.0 * christoffel[i] * velocity;
        for (size_t j = 0; j < 4; j++) {
            jac(i, 4 + j) = tmp[j];
        }
    }
    if (!proper_time) {
        double t0 = (velocity.transpose() * christoffel[0] * velocity);
        t0 /= c_c;
        for (size_t i = 0; i < 4; i++) {
            auto tmp = 2.0 * christoffel[0] * velocity * (velocity[i] / c_c);
            for (size_t j = 0; j < 4; j++) {
                jac(i, 4 + j) += tmp[j];
            }

            jac(i, i) += t0;
        }
    }
    return std::make_tuple(acc4, jac);


}