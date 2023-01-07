#pragma once
#include <tuple>
#include <Eigen/Eigen>
#include <vector>
using namespace Eigen;
using namespace std;



class Newton {

public:
	Newton() { clear(); }
	void clear() { this->masses.clear(); }
	void add_mass(double mass, Vector3d r) { tuple<double,Vector3d> gg = make_tuple(mass, r); this->masses.push_back(gg);}



	std::tuple<double, int> step_implicit_midpoint(Vector3d& x, Vector3d& v, const std::tuple<Vector3d, Vector3d>& x_v_old, const double& dt, const double tol = 1e-8);
	void step_rk4(Vector3d& x, Vector3d& v, const std::tuple<Vector3d, Vector3d>& x_v_old, const double& dt);

	Vector3d calculate_acc_newton(const Vector3d& x);
	tuple<Vector3d, Matrix3d> calculate_acc_jac_newton(const Vector3d& x);
	tuple<Vector3d, Matrix3d> calculate_acc_jac_newton_numerical(const Vector3d& xx,double eps=1e-2);

private:
	void implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& now, const VectorXd& old, const double& dt);
	
	Vector<double, 6> rk_f(const Vector<double, 6>& d);

	vector<tuple<double, Vector3d>> masses;
};

Vector3d Newton::calculate_acc_newton(const Vector3d& xx) {
	double g = 6.67430e-11;

	Vector3d acc = { 0, 0, 0 };
	for (auto it = this->masses.begin(); it != this->masses.end(); it++) {
		Vector3d r0 = get<1>(*it);
		double mass = get<0>(*it);

		Vector3d r = xx - r0;
		double r2 = r.transpose().dot(r);
		Vector3d er = r / sqrt(r2);
		acc +=  -mass * g * er / r2;
	}
	return acc;
}


void Newton::step_rk4(Vector3d& x, Vector3d& v, const std::tuple<Vector3d, Vector3d>& x_v_old, const double& dt) {
	Vector3d x_old = std::get<0>(x_v_old);
	Vector3d v_old = std::get<1>(x_v_old);
	Vector<double, 6> old, newd;
	for (size_t i = 0; i < 3; i++) {
		old[i] = x_old[i];
		old[i + 3] = v_old[i];
	}

	Vector<double, 6> k1, k2, k3, k4;
	k1 = rk_f(old);
	k2 = rk_f(old + 0.5 * dt * k1);
	k3 = rk_f(old + 0.5 * dt * k2);
	k4 = rk_f(old + dt * k3);

	newd = old + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

	for (size_t i = 0; i < 3; i++) {
		x[i] = newd[i];
		v[i] = newd[3 + i];
	}

}

Vector<double, 6> Newton::rk_f(const Vector<double, 6>& d) {
	Vector<double, 6> tmp;
	Vector<double, 3> x = { d[0],d[1],d[2]};
	Vector<double, 3> v = { d[3],d[4],d[5]};

	tmp << v, calculate_acc_newton(x);
	return tmp;
}

tuple<Vector3d, Matrix3d> Newton::calculate_acc_jac_newton_numerical(const Vector3d& xx,double eps) {
	Vector3d acc = calculate_acc_newton(xx);
	Matrix3d jac = Matrix3d::Zero();

	for (size_t i = 0; i < 3; i++) {
		Vector3d xp = 1.0 * xx;
		Vector3d xm = 1.0 * xx;
		xp[i] += eps;
		xm[i] -= eps;

		Vector3d dacc = 0.5 * (calculate_acc_newton(xp) - calculate_acc_newton(xm)) / eps;
		jac(0, i) = dacc(0);
		jac(1, i) = dacc(1);
		jac(2, i) = dacc(2);
	}
	return make_tuple(acc, jac);
}
tuple<Vector3d,Matrix3d> Newton::calculate_acc_jac_newton(const Vector3d& xx) {
	double g = 6.67430e-11;

	Vector3d acc = { 0, 0, 0 };
	Matrix3d accjac = Matrix3d::Zero();
	for (auto it = this->masses.begin(); it != this->masses.end(); it++) {
		Vector3d r0 = get<1>(*it);
		double mass = get<0>(*it);

		Vector3d r = xx - r0;
		double r2 = r.transpose().dot(r);
		double rr = sqrt(r2);
		double r3 = r2 * rr;
		double r5 = r2 * r3;
		Vector3d er = r / rr;
		acc += -mass * g * er / r2;


		accjac(0, 0) += -g * mass / r3;
		accjac(0, 0) += 1.5 * g * mass * r(0) * 2 * r(0)/r5;
		accjac(0, 1) += 1.5 * g * mass * r(0) * 2 * r(1) / r5;
		accjac(0, 2) += 1.5 * g * mass * r(0) * 2 * r(2) / r5;

		accjac(1, 1) += -g * mass / r3;
		accjac(1, 1) += 1.5 * g * mass * r(1) * 2 * r(1) / r5;
		accjac(1, 0) += 1.5 * g * mass * r(1) * 2 * r(0) / r5;
		accjac(1, 2) += 1.5 * g * mass * r(1) * 2 * r(2) / r5;

		accjac(2, 2) += -g * mass / r3;
		accjac(2, 2) += 1.5 * g * mass * r(2) * 2 * r(2) / r5;
		accjac(2, 0) += 1.5 * g * mass * r(2) * 2 * r(0) / r5;
		accjac(2, 1) += 1.5 * g * mass * r(2) * 2 * r(1) / r5;

	}
	return make_tuple(acc,accjac);
}
void Newton::implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& nowd, const VectorXd& old, const double& dt) {
	
	Vector3d xx;
	xx << 0.5 * (nowd[0] + old[0]), 0.5 * (nowd[1] + old[1]), 0.5 * (nowd[2] + old[2]);

	auto [acc, accjac] = calculate_acc_jac_newton(xx);
	for (size_t i = 0; i < 3; i++) {
		res[i] = (nowd[i] - old[i]) / dt - (nowd[i + 3] + old[i + 3]) / 2;
		res[i + 3] = (nowd[i + 3] - old[i + 3]) / dt - acc[i];
	}
	for (size_t i = 0; i < 3; i++) {
		J(i, i) = 1.0/dt;
		J(i, i + 3) = -0.5;

		J(i + 3, i + 3) = 1.0 / dt;
		J(i + 3, 0) = -0.5 * accjac(i, 0);
		J(i + 3, 1) = -0.5 * accjac(i, 1);
		J(i + 3, 2) = -0.5 * accjac(i, 2);
	}
	
}


std::tuple<double, int> Newton::step_implicit_midpoint(Vector3d& x, Vector3d& v, const std::tuple<Vector3d, Vector3d>& x_v_old, const double& dt, const double tol) {

	VectorXd dx, errs;

	VectorXd resids;
	MatrixXd J;

	resids = Vector<double, 6>::Zero();
	J = Matrix<double, 6, 6>::Zero();

	auto [x_old, v_old] = x_v_old;
	Vector<double, 6> nowd, oldd;
	for (size_t i = 0; i < 3; i++) {
		oldd[i] = x_old[i];
		oldd[i + 3] = v_old[i];
	}


	//
	// Init guess - Semi implicit Euler
	////


	Vector3d acc = calculate_acc_newton(x_old);
	for (size_t i = 0; i < 3; i++) {
		nowd[3 + i] = v_old[i] + dt * acc[i];
		nowd[i] = x_old[i] + dt * nowd[3 + i];
	}



	double err;
	//Newtons method
	int n = 1;
	for (size_t i = 0; i < 10; i++) {


		implicit_midpoint_resids(resids, J, nowd, oldd, dt);

		dx = J.llt().solve(-resids);
		nowd = nowd + dx;

		errs = J * dx - resids;
		err = -1.0;
		for (size_t j = 0; j < 6; j++) {
			if (err < abs(errs[j]))
				err = abs(errs[j]);

		}

		if (err < tol)
			break;
		n += 1;
	}
	for (size_t i = 0; i < 3; i++) {
		x[i] = nowd[i];
		v[i] = nowd[3 + i];
	}
	return std::make_tuple(err, n);
};

