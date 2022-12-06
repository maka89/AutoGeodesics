#pragma once
#include "autogeodesics.h"
void AutoGeodesics::calculate_acc_velocity_jacobian(Vector4d& acc, Matrix4d& J, const Vector4d& velocity, const Vector4d& x)
{


std::array<Matrix4d, 4> christoffel = calculate_christoffel(x);
acc = get_acc(christoffel, velocity);

for (size_t i = 0; i < 4; i++)
	J.row(i) = 2.0 * christoffel[i] * velocity;

}

int AutoGeodesics::velocity_verlet_resids(Vector4d &resids, Matrix4d &jacc, const Vector4d& velocity, const std::array<Matrix4d, 4>& chr, const Vector4d& v_old, const Vector4d& acc_old, const Vector4d& x, const double& dt) {

	resids = (velocity - v_old) - 0.5 * (acc_old + get_acc(chr, velocity)) * dt;

	for (size_t i = 0; i < 4; i++)
		jacc.row(i) = chr[i] * velocity * dt;
	jacc = jacc+ Matrix4d::Identity();

	return 0;
}

void AutoGeodesics::step_euler(Vector4d& x, Vector4d& v, Vector4d& acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old,  const double& dt) {
	auto [x_old, v_old, acc_old] = x_v_a_old;
	v = v_old + acc_old * dt;
	x = x_old + v_old * dt;
	acc = calculate_acc(x, v);
}
void AutoGeodesics::step_euler_cromer(Vector4d& x, Vector4d& v, Vector4d& acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old,  const double& dt) {
	auto [x_old, v_old, acc_old] = x_v_a_old;
	v = v_old + acc_old * dt;
	x = x_old + v * dt;
	acc = calculate_acc(x, v);
}
void AutoGeodesics::step_euler_midpoint(Vector4d& x, Vector4d& v, Vector4d& acc, const std::tuple<Vector4d, Vector4d, Vector4d>& x_v_a_old, const double& dt) {
	auto [x_old, v_old, acc_old] = x_v_a_old;
	v = v_old + acc_old * dt;
	x = x_old + 0.5*(v+v_old) * dt;
	acc = calculate_acc(x, v);
}

double AutoGeodesics::step_velocity_verlet(Vector4d &x,Vector4d &v,Vector4d &acc, const std::tuple<Vector4d, Vector4d, Vector4d> &x_v_a_old,  const double& dt,const double tol) {

	Vector4d dv, errs, resids;
	Matrix4d J;

	auto [x_old, v_old, acc_old] = x_v_a_old;

	//Update position
	x = x_old + v_old * dt + 0.5 * acc_old * dt * dt;

	
	v = v_old +  acc_old* dt; //init guess for velocity;


	std::array<Matrix4d, 4> christoffel = calculate_christoffel(x);

	//Newtons method
	double err = 0.0;
	for (size_t i = 0; i < 30; i++) {



		velocity_verlet_resids(resids, J, v, christoffel, v_old, acc_old, x, dt);

		dv = J.llt().solve(-resids);
		v = v + dv;

		errs = J * dv - resids;
		err = -1.0;
		for (size_t j = 0; j < 4; j++) {
			if (err < abs(errs[j])) 
				err = abs(errs[j]);
			
		}

		if (err < tol) 
			break;
		
	}

	//Update velocity and acceleration
	acc = get_acc(christoffel, v);


	return err;
};

