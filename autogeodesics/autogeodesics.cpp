#pragma once
#include "autogeodesics.h"
#include <iostream>
void AutoGeodesics::implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& nowd, const VectorXd& old, const double& dt) {

	Vector<double, 8> now;
	for (size_t i = 0; i < 8; i++)
		now[i] = nowd[i];
	Vector4d xnow = { now[0],now[1],now[2],now[3] };
	Vector4d xold = { old[0],old[1],old[2],old[3] };
	Vector4d acc;

	Vector4d u = { 0.5 * (now[4] + old[4]),0.5 * (now[5] + old[5]),0.5 * (now[6] + old[6]),0.5 * (now[7] + old[7]) };

	Vector<double, 8> resids;
	if (this->proper_time) {
		for (size_t i = 0; i < 4; i++)
			resids[i] = now[i] - old[i] - 0.5 * dt * (now[i + 4] + old[i + 4]);

		auto [acc, jacc] = calculate_acc_jac(0.5 * (now + old));


		for (size_t i = 0; i < 4; i++)
			resids[4 + i] = now[4 + i] - old[4 + i] - acc[i] * dt;

		for (size_t i = 0; i < 4; i++) {
			J(i, i) = 1.0;
			J(i, i + 4) = -0.5 * dt;
		}
		for (size_t i = 0; i < 4; i++) {
			J.row(4 + i) = -jacc.row(i)*dt;
			J(4 + i, 4 + i) += 1.0;
		}

	}
	else {


		for (size_t i = 0; i < 4; i++)
			resids[i] = now[i] - old[i] - 0.5 * dt * (now[i + 4] + old[i + 4]);

		auto[acc,jacc] = calculate_acc_jac(0.5 * (now+old));

		for (size_t i = 0; i < 4; i++)
			resids[4 + i] = now[4 + i] - old[4 + i] - acc[i] * dt;
		
		J = Matrix<double, 8, 8>::Zero();
		for (size_t i = 0; i < 4; i++) {
			J(i, i) = 1.0;
			J(i, i + 4) = -0.5 * dt;
		}
		for (size_t i = 0; i < 4; i++) {
			J.row(4 + i) = -jacc.row(i)*dt;
			J(4 + i, 4 + i) += 1.0;
		}
	}
	res = resids;
}

Vector<double,8> AutoGeodesics::rk_f(const Vector<double, 8> &d) {
	Vector<double, 8>tmp;
	Vector4d x = { d[0],d[1],d[2],d[3] };
	Vector4d v = { d[4],d[5],d[6],d[7] };

	tmp(seq(0, 4)) = v;
	tmp(seq(4, last)) = calculate_acc(x, v);
	return tmp;
}
void AutoGeodesics::step_rk4(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt) {
	Vector4d x_old = std::get<0>(x_v_old);
	Vector4d v_old = std::get<1>(x_v_old);
	Vector<double, 8> old,newd;
	for (size_t i = 0; i < 4; i++) {
		old[i] = x_old[i];
		old[i + 4] = v_old[i];
	}

	Vector<double, 8> k1, k2, k3, k4;
	k1 = rk_f(old);
	k2 = rk_f(old + 0.5 * dt * k1);
	k3 = rk_f(old + 0.5 * dt * k2);
	k4 = rk_f(old + dt * k3);

	newd = old+dt*(k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

	for (size_t i = 0; i < 4; i++) {
		x[i] = newd[i];
		v[i] = newd[4 + i];
	}



}

std::tuple<double,int> AutoGeodesics::step_implicit_midpoint(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt, const double tol) {

	VectorXd dx, errs;

	VectorXd resids,resids2;
	MatrixXd J,J2;

	resids = Vector<double, 8>::Zero();
	J = Matrix<double, 8, 8>::Zero();

	resids2 = Vector<double, 8>::Zero();
	J2 = Matrix<double, 8, 8>::Zero();

	Vector4d x_old = std::get<0>(x_v_old);
	Vector4d v_old = std::get<1>(x_v_old);

	Vector<double, 8> nowd, oldd;
	for (size_t i = 0; i < 4; i++){
		oldd[i] = x_old[i];
		oldd[i + 4] = v_old[i];
	}


	//
	// Init guess - Semi implicit Euler
	////
	
	Vector4d acc = calculate_acc(x_old,v_old);
	for (size_t i = 0; i < 4; i++) {
		nowd[4 + i] = v_old[i] + dt * acc[i];
		nowd[i] = x_old[i]+ dt * nowd[4 + i];
	}
	

	
	double err;
	//Newtons method
	int n = 1;
	for (size_t i = 0; i < 10; i++) {

		err = -1.0;
		for (size_t j = 0; j < 8; j++) {
			if (err < abs(resids[j]))
				err = abs(resids[j]);

		}
		if (err < tol && i >0)
			break;

		implicit_midpoint_resids(resids, J, nowd, oldd, dt);

		dx = -J.inverse() * resids;

		nowd = nowd + dx;

		

		
		n += 1;
	}
	for (size_t i = 0; i < 4; i++) {
		x[i] = nowd[i];
		v[i] = nowd[4 + i];
	}
	return std::make_tuple(err,n);
};


/*

void AutoGeodesics::implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& nowd, const VectorXd& old, const double& dt) {

	Vector<var, 8> now;
	for (size_t i = 0; i < 8; i++)
		now[i] = nowd[i];
	Vector4var xnow = { now[0],now[1],now[2],now[3] };
	Vector4var xold = { old[0],old[1],old[2],old[3] };
	Vector4var acc;

	Vector4var u = { 0.5 * (now[4] + old[4]),0.5 * (now[5] + old[5]),0.5 * (now[6] + old[6]),0.5 * (now[7] + old[7]) };

	Vector<var, 8> resids;
	if (this->proper_time) {
		for (size_t i = 0; i < 4; i++)
			resids[i] = now[i] - old[i] - 0.5 * dt * (now[i + 4] + old[i + 4]);

		acc = calculate_acc(0.5 * (xnow + xold), u);
		for (size_t i = 0; i < 4; i++)
			resids[4 + i] = now[4 + i] - old[4 + i] - acc[i] * dt;

	}
	else {


		resids[0] = now[0] - old[0] - c_c * dt;
		for (size_t i = 1; i < 4; i++)
			resids[i] = now[i] - old[i] - 0.5 * dt * (now[i + 4] + old[i + 4]);

		acc = calculate_acc(0.5 * (xnow + xold), u);

		resids[4] = now[4] - old[4];
		for (size_t i = 1; i < 4; i++)
			resids[4 + i] = now[4 + i] - old[4 + i] - acc[i] * dt;

	}
	res = resids.cast<double>();
	for (size_t i = 0; i < 8; i++)
		J.row(i) = gradient(resids[i], now);
}
*/