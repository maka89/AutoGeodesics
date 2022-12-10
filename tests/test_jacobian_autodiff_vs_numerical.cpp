#include <autogeodesics.h>
#include <iostream>

Matrix4var Metricfn(const Vector4var& x) {
	return AutoGeodesics::Metrics::schwarzschild_cartesian(x, 5.97e24, Vector3d({ 0.0,0.0,0.0 }));
}

int compareJ(Matrix<double, 4, 8> J_num, Matrix<double, 4, 8> J_auto) {
	for (size_t i = 0; i < 4; i++) {
		std::cout << "X[" << i << "]" << std::endl;
		std::cout << "Num:  " << J_num.row(i) << std::endl;
		std::cout << "Auto: " << J_auto.row(i) << std::endl << std::endl;
	}
	return 0;
}



int test(Vector3d x0,Vector3d vel3,Vector<double,8> dx) {


	AutoGeodesics ag = AutoGeodesics(false);
	ag.setMetFn(&Metricfn);

	Vector4d x = Vector4d({ 0.0,x0[0],x0[1],x0[2] });
	Vector4d vel = ag.setup_fourvelocity(x,vel3);
	
	Vector<double, 8> u;
	u << x, vel;
	
	auto t = ag.calculate_acc_jac(u);
	Matrix<double, 4, 8> J = std::get<1>(t);
	Matrix<double, 4, 8> J2 = Matrix<double, 4, 8>();

	for (size_t j = 0; j < 8; j++) {
		Vector<double,8> up = 1.0 * u;
		up[j] += dx[j];
		Vector<double, 8> um = 1.0 * u;
		um[j] -= dx[j];

		auto [acp, _p] = ag.calculate_acc_jac(up);
		auto [acm, _m] = ag.calculate_acc_jac(um);

		J2.col(j) = 0.5 * (acp - acm) / dx[j];

	}
	compareJ(J2, J);
	return 0;
}

int main() {
	double alpha = 1e-23;
	double beta = 1e2;
	Vector3d v,x;
	Vector<double, 8> dx;

	std::cout << "X= 6371km from earth" << std::endl;
	x = { 6371000.0,0.0,0.0 };
	v = { 0.0,100000.0,0.0 };
	dx = { 1e-1, 1e-4, 1e-4, 1e-4, 1e-1, beta + alpha * v[0], beta + alpha * v[1], beta+ alpha * v[2] };
	test(x,v, dx); 



	std::cout << "X= 1000.0 km from earth" << std::endl;
	x = { 0.02,0.02,0.0 };
	v = { 10000.0,10000.0,0.0 };
	dx = { 1e-1, 1e-4, 1e-4, 1e-4, 1e-1, beta + alpha * v[0], beta + alpha * v[1], beta + alpha * v[2] };
	test(x, v, dx);
}