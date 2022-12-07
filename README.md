# AutoGeodesics
Easily integrate the geodesics equation using automatic differentiation. C++17

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

### Set a metric function
Write your own metric function or use built in templates. 
 ~~~c++
 
 //Write your own...
Matrix4real schwarzschild(const Vector4real& x) {
    double mass = 5.972e24;
    double rs = 2.0 * c_g * mass / pow(c_c, 2);
    Matrix4real metric = Matrix4real::Zero();
    metric(0, 0) = -(1.0 - rs / x[1]);
    metric(1, 1) = 1.0 / (1.0 - rs / x[1]);
    metric(2, 2) = x[1] * x[1];
    metric(3, 3) = x[1] * x[1] * sin(x[2]) * sin(x[2]);

    return -1 * metric;
}

//Or use built-in templates
Matrix4real schwarzschild_cartesian(const Vector4real& x) {
    double mass = 5.972e24;
    Vector3d x0 = { 0.0,0.0,0.0 };
    return AutoGeodesics::Metrics::schwarzchild_cartesian(x, mass, x0);
}

AutoGeodesics ag = AutoGeodesics();
ag.setMetFn(&schwarzschild_cartesian);
~~~

### Calculate acceleration
Calculate the 4-acceleration for a set of velocity and position.
 ~~~c++


Vector4d x = {0.0,6371000.0,0.0,0.0};
Vector3d vel3 = {0.0,sqrt(9.81*6371000.0),0.0}

Vector4d velocity = ag.setup_fourvelocity(x, vel3); //Turn velocity into 4-velocity.
Vector4d acc = ag.calculate_acc(x, velocity);
~~~  
### Integrate the Geodesics Equation...
Built-in methods lets you calculate the next position/velocity. The Geodesic Equation can be integrated using velocity-verlet integration.

 ~~~c++
 
 int steps = 200000;
 double t_end = 2*3.1415*sqrt(6371000.0/9.81);
 double dt = t_end / steps;

Vector4d xo, velo, acco; //Last timestep position/velocity/acceleration.
    
 for (int i = 0; i < steps; i++) {
     acco = acc;
     velo = vel;
     xo = x;
     ag.step_velocity_verlet(x, vel, acc, std::make_tuple(xo, velo, acco), dt, 1e-6); //x,vel,acc is passed by reference and overwritten.
 }
 ~~~   

### Use your own numerical methods
Common integration schemes , like leapfrog and velocity verlet, require implicit methods for calculating the velocity  $v_{i+1} = f(v_{i+1})$. 
Efficiently solving these equations, require the acceleration and the jacobian of the acceleration wrt. the velocity.
The class has methods for calculating the jacobian.
 ~~~c++
Vector4d acc;
Matrix4d J;
ag.calculate_acc_velocity_jacobian(acc, J, velocity, x); //acc,J passed by reference and overwritten.
~~~
