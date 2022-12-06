#include <autogeodesics.h>
#include <iostream>

using namespace std;

Matrix4real schwarzchild_cartesian(const Vector4real& x) {
	double mass = 5.972e24;
	Vector3d x0 = { 0.0,0.0,0.0 };
	return AutoGeodesics::Metrics::schwarzchild_cartesian(x, mass, x0);
}

int main() {
	AutoGeodesics ag = AutoGeodesics();
	ag.setMetFn(&schwarzchild_cartesian);
	Vector4d x = { 0.0,6371000.0,0.0,0.0 };

	Vector4d velocity = ag.setup_fourvelocity(x, Vector3d({ 0.0,sqrt(9.81997 * 6371000.0),0.0 })); //Turn velocity into 4-velocity.
	Vector4d acc = ag.calculate_acc(x, velocity);

	cout << acc << endl;

	int steps = 200;
	double t_end = 2 * 3.14159265 * sqrt(6371000.0 / 9.81997);
	double dt = t_end / steps;

	Vector4d xo, velo, acco; //Last timestep position/velocity/acceleration.


	cout << x[0] / c_c << ", " << x[1] << ", " << x[2] << ", " << x[3] << endl;
	for (int i = 0; i < steps; i++) {
		acco = acc;
		velo = velocity;
		xo = x;
		ag.step_velocity_verlet(x, velocity, acc, std::make_tuple(xo, velo, acco), dt, 1e-8); //x,vel,acc is passed by reference and overwritten.

		cout << x[0] / c_c << ", " << x[1] << ", " << x[2] << ", " << x[3] << endl;
	}


}