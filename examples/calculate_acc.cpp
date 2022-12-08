// C++ includes
#include <iostream>
#include <autogeodesics.h>


inline Matrix4var schwarzschild_cart(const Vector4var& x) {
    double mass = 5.972e24;
    Vector3d x0 = Vector3d::Zero();

    return AutoGeodesics::Metrics::schwarzschild_cartesian(x, mass, x0);
}


int main() {
    Vector4d x,velocity;

    // Setup. proper_time=false. Use coordinate time.
    AutoGeodesics ag = AutoGeodesics(false);

    //Set metric function.
    ag.setMetFn(schwarzschild_cart); 

    // Position y=Radius of Earth.
    x = { 0.0,0.0,6371000.0,0.0 }; 

    // Function to setup all 4 components of velocity. Mostly useful for proper_time=true.
    // In this case velocity could be set up with velocity << c,0,0,0;
    velocity = ag.setup_fourvelocity(x, Vector3d({ 0.0,0.0,0.0 })); 


    //Calculates acceleration. (0 component included).
    std::cout << ag.calculate_acc(x, velocity).cast<double>() << std::endl; 

}
