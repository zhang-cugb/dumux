## Problem setups and parameter distributions

In the main file we have solved a single-phase and a tracer simulation, which are
based on different available models in Dumux. To set up a simulation in Dumux, you have
to define compile-time `properties` and you have to implement a `Problem` class containing
initial and boundary conditions, as well as a `SpatialParameters` class defining the distributions
of porous medium parameters such as e.g. permeability and porosity. Please see the documentations
provided in the following to see how this is realized in this example:

### Single-phase flow simulation

* Property definitions: __TODO: LINK DOCU OF `properties_1p.hh`__
* Problem class: __TODO: LINK DOCU OF `problem_1p.hh`__
* Spatial parameters class: __TODO: LINK DOCU OF `spatialparams_1p.hh`__

### Tracer transport simulation

* Property definitions: __TODO: LINK DOCU OF `properties_tracer.hh`__
* Problem class: __TODO: LINK DOCU OF `problem_tracer.hh`__
* Spatial parameters class: __TODO: LINK DOCU OF `spatialparams_tracer.hh`__
