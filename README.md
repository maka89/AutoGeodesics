# AutoGeodesics - API2

This branch does the same as the main branch, but uses forward mode autodiff instead of reverse mode.
This gives a 5x speedup for implicit midpoint method.
A separate function has to be specified for each of the 10 unique components of the metric.

~~~c++
AutoGeodesics ag = AutoGeodesics(false);
ag.setMetFnComp(schwarzschild_cartesian_00, 0);

//Zero components needn't be specified.
//ag.setMetFnComp(schwarzschild_cartesian_01, 1);
//ag.setMetFnComp(schwarzschild_cartesian_02, 2);
//ag.setMetFnComp(schwarzschild_cartesian_03, 3);

ag.setMetFnComp(schwarzschild_cartesian_11, 4);
ag.setMetFnComp(schwarzschild_cartesian_12, 5);
ag.setMetFnComp(schwarzschild_cartesian_13, 6);
ag.setMetFnComp(schwarzschild_cartesian_22, 7);
ag.setMetFnComp(schwarzschild_cartesian_23, 8);
ag.setMetFnComp(schwarzschild_cartesian_33, 9);
~~~

See tests/examples for more info on how to use.

