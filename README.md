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
Matrix4var schwarzschild(const Vector4var& x) {
    double mass = 5.972e24;
    double rs = 2.0 * c_g * mass / pow(c_c, 2);
    Matrix4var metric = Matrix4var::Zero();
    metric(0, 0) = -(1.0 - rs / x[1]);
    metric(1, 1) = 1.0 / (1.0 - rs / x[1]);
    metric(2, 2) = x[1] * x[1];
    metric(3, 3) = x[1] * x[1] * sin(x[2]) * sin(x[2]);

    return -1 * metric;
}

//Or use built-in templates
Matrix4var schwarzschild_cart(const Vector4var& x) {
    double mass = 5.972e24;
    Vector3d x0 = Vector3d::Zero();

    return AutoGeodesics::Metrics::schwarzschild_cartesian(x, mass, x0);
}

int main(){

Vector4d x,velocity;

// Setup. proper_time=false. Use coordinate time.
AutoGeodesics ag = AutoGeodesics(false);

//Set metric function.
ag.setMetFn(schwarzschild_cart); 
~~~

### Calculate acceleration
Calculate the 4-acceleration for a set of velocity and position.
 ~~~c++


 // Position y=Radius of Earth.
 x = { 0.0,0.0,6371000.0,0.0 }; 

 // Function to setup all 4 components of velocity. Mostly useful for proper_time=true.
 velocity = ag.setup_fourvelocity(x, Vector3d({ 0.0,0.0,0.0 })); 

 //Calculates acceleration. (0 component included).
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
