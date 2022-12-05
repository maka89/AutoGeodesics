# AutoRelativity
Easily integrate the geodesics equation using automatic differentiation.

## Introduction
The purpose of this project is to make it easy to calculate paths of test particles being affected by gravity under the theory of general relativity. Using automatic differentiation, the user can specify a metric, without having to worry about its derivatives, which are needed to calculate the geodesics.

The project aims to do these things:
- Calculates the acceleration of the test particle, given a metric and initial position/velocity.
- Automatically calculates the next velocity/position for some popular integration schemes, such as velocity-verlet.

## Dependencies
Two header-only libraries:
- https://github.com/autodiff/autodiff
- Eigen
## Usage
### Calculate acceleration
Calculate the 4-acceleration for a metric and initial velocity/position.
 ~~~c++
//Write your own metric function, or use built-in templates.
Matrix4real schwarzchild_cartesian(const Vector4real& x) {
    double mass = 5.972e24;
    Vector3d x0 = { 0.0,0.0,0.0 };
    return AutoGeodesic::Metrics::schwarzchild_cartesian(x, mass, x0);
}

AutoGeodesic ag = AutoGeodesic();
ag.setMetFn(&schwarzchild_cartesian);

Vector4d x = {0.0,6371000.0,0.0,0.0};
Vector3d vel3 = {0.0,sqrt(9.81*6371000.0),0.0}

Vector4d vel = ag.setupfourvelocity(x,vel3); //Turn velocity into 4-velocity.
Vector4d acc = ag.calculate_acc(x,vel);
~~~  
### Integrate the Geodesics Equation...
Built-in methods lets you calculate the next position/velocity. The Geodesic Equation can be integrated using velocity-verlet integration.

 ~~~c++
 int steps = 200000;
 double t_end = 2*3.1415*sqrt(6371000.0)/9.81;
 double dt = t_end / steps;

Vector4d x, vel, acc, xo, velo, acco;
    
 for (int i = 0; i < steps; i++) {
     acco = acc;
     velo = vel;
     xo = x;
     ag.step_velocity_verlet(x, vel, acc, std::make_tuple(xo, velo, acco), dt, 1e-6);
 }
 ~~~   
