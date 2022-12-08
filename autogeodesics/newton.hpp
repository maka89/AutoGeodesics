#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include <tuple>
#include <Eigen/Eigen>
using namespace autodiff;
using namespace Eigen;




class Newton {

public:
    Newton() :potential(NULL) {};
   
	Vector4d calculate_acceleration_newton(const Vector4d &x);
	void set_potential(var(*pot)(const Vector4var& x)) { this->potential = pot; }


	std::tuple<double, int> step_implicit_midpoint(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt, const double tol = 1e-8);
	void step_rk4(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt);

private:
	void implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& now, const VectorXd& old, const double& dt);
	Vector4var calculate_acc_newton(const Vector4var& x);
	Vector<double, 8> rk_f(const Vector<double, 8>& d);
	var(*potential)(const Vector4var& x);
};

Vector4var Newton::calculate_acc_newton(const Vector4var& xx) {
	double g = 6.67430e-11;

	var potential = -g*this->potential(xx);
	Vector4var acc;
	auto [a,b,c,d] = derivatives(potential, wrt(xx[0],xx[1],xx[2],xx[3]));
	acc << a, b, c, d;
	

    return acc;
}

Vector4d Newton::calculate_acceleration_newton(const Vector4d &x) {
	Vector4var vv;
	vv << x[0], x[1], x[2], x[3];
	return calculate_acc_newton(vv).cast<double>();
}

void Newton::step_rk4(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt) {
	Vector4d x_old = std::get<0>(x_v_old);
	Vector4d v_old = std::get<1>(x_v_old);
	Vector<double, 8> old, newd;
	for (size_t i = 0; i < 4; i++) {
		old[i] = x_old[i];
		old[i + 4] = v_old[i];
	}

	Vector<double, 8> k1, k2, k3, k4;
	k1 = rk_f(old);
	k2 = rk_f(old + 0.5 * dt * k1);
	k3 = rk_f(old + 0.5 * dt * k2);
	k4 = rk_f(old + dt * k3);

	newd = old + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

	for (size_t i = 0; i < 4; i++) {
		x[i] = newd[i];
		v[i] = newd[4 + i];
	}

}

Vector<double, 8> Newton::rk_f(const Vector<double, 8>& d) {
	Vector<double, 8>tmp;
	Vector4d x = { d[0],d[1],d[2],d[3] };
	Vector4d v = { d[4],d[5],d[6],d[7] };

	tmp(seq(0, 4)) = v;
	tmp(seq(4, last)) = calculate_acceleration_newton(x);
	return tmp;
}

void Newton::implicit_midpoint_resids(VectorXd& res, MatrixXd& J, const VectorXd& nowd, const VectorXd& old, const double& dt) {

	Vector<var, 8> now;
	for (size_t i = 0; i < 8; i++)
		now[i] = nowd[i];
	Vector4var xnow = { now[0],now[1],now[2],now[3] };
	Vector4var xold = { old[0],old[1],old[2],old[3] };
	Vector4var acc;

	Vector4var u = { 0.5 * (now[4] + old[4]),0.5 * (now[5] + old[5]),0.5 * (now[6] + old[6]),0.5 * (now[7] + old[7]) };

	Vector<var, 8> resids;

	resids[0] = now[0] - old[0];
	for (size_t i = 1; i < 4; i++)
		resids[i] = now[i] - old[i] - 0.5 * dt * (now[i + 4] + old[i + 4]);

	acc = calculate_acc_newton(0.5 * (xnow + xold));

	resids[4] = now[4] - old[4];
	for (size_t i = 1; i < 4; i++)
		resids[4 + i] = now[4 + i] - old[4 + i] - acc[i] * dt;

	
	res = resids.cast<double>();
	for (size_t i = 0; i < 8; i++)
		J.row(i) = gradient(resids[i], now);
}


std::tuple<double, int> Newton::step_implicit_midpoint(Vector4d& x, Vector4d& v, const std::tuple<Vector4d, Vector4d>& x_v_old, const double& dt, const double tol) {

	VectorXd dx, errs;

	VectorXd resids;
	MatrixXd J;

	resids = Vector<double, 8>::Zero();
	J = Matrix<double, 8, 8>::Zero();

	auto [x_old, v_old] = x_v_old;
	Vector<double, 8> nowd, oldd;
	for (size_t i = 0; i < 4; i++) {
		oldd[i] = x_old[i];
		oldd[i + 4] = v_old[i];
	}


	//
	// Init guess - Semi implicit Euler
	////


	Vector4d acc = calculate_acc_newton(x_old.cast<var>()).cast<double>();
	for (size_t i = 0; i < 4; i++) {
		nowd[4 + i] = v_old[i] + dt *acc[i] ;
		nowd[i] = x_old[i] + dt * nowd[4 + i];
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
		for (size_t j = 0; j < 8; j++) {
			if (err < abs(errs[j]))
				err = abs(errs[j]);

		}

		if (err < tol)
			break;
		n += 1;
	}
	for (size_t i = 0; i < 4; i++) {
		x[i] = nowd[i];
		v[i] = nowd[4 + i];
	}
	return std::make_tuple(err, n);
};

