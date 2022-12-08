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
//Calculate 1 circular orbit at earth radius. Use 1000 steps.
int steps = 1000;
double dt= 2 * 3.141592 * sqrt(x[2] / 9.81997) / steps;
Vector4d velo,xo;
for (int i = 0; i < steps; i++) {
    velo = velocity;
    xo = x;


    // Make one step using RungeKutta4
    ag.step_rk4(x, velocity, std::make_tuple(xo, velo), dt);

    //Alternatively, use implicit midpoint rule
    //auto [err, niter] = ag.step_implicit_midpoint(x, velocity, std::make_tuple(xo, velo), dt, 1e-3);


    cout << x[0]/c_c<<", "<<x[1]<<", "<<x[2]<<", "<<x[3] << endl;

}
 ~~~   

### Use your own numerical methods
Common integration schemes , like leapfrog and velocity verlet, require implicit methods for calculating the velocity  $v_{i+1} = f(v_{i+1})$. 
Efficiently solving these equations, require the acceleration and the jacobian of the acceleration wrt. the velocity.
The class has methods for calculating the jacobian.

The function `Vector4var AutoGeodesics::calculate_acc(const Vector4var& x, const Vector4var& velocity)` function uses c++ autodiff functions. The Jacobian of the acceleration can be calculated wrt. positiono or velocity

~~~

J=Matrix4d;
Vector4var acc = calculate_acc(x,v);

//Jacobian (acc wrt x)
for (size_t i = 0; i < 4; i++)
  J.row(i) = gradient(resids[i], x);
  
J=Matrix4d;
// Jacobian(acc wrt v)
for (size_t i = 0; i < 4; i++)
		J.row(i) = gradient(resids[i], v);
  
Vector4var acc = calculate_acc(const Vector<var,8> x){
  return calculate_acc(x(seq(0,4),x(seq(4,last));
}

J=Matrix<double,4,8>
//Jacobian (acc wrt x)

Vector<var,8> xx;
xx << x[0],x[1],x[2],x[3],v[0],v[1],v[2],v[3];
for (size_t i = 0; i < 4; i++)
  J.row(i) = gradient(resids[i], xx);
~~~

