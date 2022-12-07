#include <autogeodesics.h>
#include <iostream>
Matrix4real schwarzschild_cartesian(const Vector4real& x) {
	double mass = 5.972e24;
	Vector3d x0 = { 0.0,0.0,0.0 };
	return AutoGeodesics::Metrics::schwarzschild_cartesian(x, mass, x0);
}

int main() {
	AutoGeodesics ag = AutoGeodesics();
	ag.setMetFn(&schwarzchild_cartesian);
	Vector4d x = { 0.0,6371000.0,0.0,0.0 };

	Vector4d velocity = ag.setup_fourvelocity(x, Vector3d({ 0.0,sqrt(9.81 * 6371000.0),0.0 })); //Turn velocity into 4-velocity.
	Vector4d acc = ag.calculate_acc(x, velocity);

	std::cout << acc << std::endl;

}
