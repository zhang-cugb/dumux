<!-- Important: This file has been automatically generated by generate_example_docs.py. Do not edit this file directly! -->

# Rotation-symmetric pressure distribution

__In this example, you will learn how to__

* solve a rotation-symmetric problem one-dimensionally
* perform a convergence test against an analytical solution
* apply the `RotationalExtrusion` filters in [ParaView](https://www.paraview.org/) for a two-dimensional visualization of the one-dimensional results


__Result__. With the `RotationalExtrusion` and the `WarpByScalar` filter in [ParaView](https://www.paraview.org/),
the pressure distribution of this example looks as shown in the following picture:

<figure>
    <center>
        <img src="img/result.png" alt="Rotation-symmetric pressure distribution" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Result TODO.</figcaption>
    </center>
</figure>


__Table of contents__. This description is structured as follows:

[[_TOC_]]


# Problem setup

We consider a single-phase problem that leads to a rotation-symmetric pressure distribution.
The following figure illustrates the setup:

<figure>
    <center>
        <img src="img/setup.svg" alt="Rotation-symmetric setup" width="60%"/>
        <figcaption> <b> Fig.1 </b> - Setup for the rotation-symmetric problem. The pressure boundary conditions are shown by the colored lines and the simulation domain is depicted in grey.</figcaption>
    </center>
</figure>

This could, for example, represent a cross section of an injection/extraction well in a homogeneous
and isotropic porous medium, where the well with radius $`r_1`$ is cut out and the
injection/extraction pressure $`p_1`$ is prescribed as a Dirichlet boundary condition. At the outer
radius $`r_2`$, we set the pressure $`p_2`$. In the polar coordinates $`r`$ and $`\varphi`$, the
solution to this problem is independent of the angular coordinate $`\varphi`$ and can be reduced to
a one-dimensional problem in the radial coordinate $`r`$. Therefore, in this example, we want to
solve the problem on a one-dimensional computational domain as illustrated by the orange line in
the above figure.

# Mathematical model

In this example we are using the single-phase model of DuMuX, which considers Darcy's law to relate
the Darcy velocity $`\textbf v`$ to gradients of the pressure $`p`$. For an isotropic medium and
neglecting gravitational forces, this can be written as:

```math
\textbf v = - \frac{k}{\mu} \text{grad} p.
```

Here, $`k`$ is the permeability of the porous medium and $`\varrho`$ and $`\mu`$ are the density
and the dynamic viscosity of the fluid. In the model, the mass balance equation for the fluid
phase is solved:

```math
\phi \frac{\partial \varrho}{\partial t} + \text{div} \left( \varrho \textbf v \right) = 0,
```

where $`\phi`$ is the porosity of the porous medium. Let us now introduce the transformation
$`(x, y)^T = \Phi ( r, \varphi )`$ from polar into cartesian coordinates (see e.g.
[wikipedia.org](https://en.wikipedia.org/wiki/Polar_coordinate_system#Converting_between_polar_and_Cartesian_coordinates)),
and denote with

```math
\tilde{p} \left( r, \varphi \right)
    = \tilde{p} \left( r \right)
    = p \left( \Phi^{-1}(x, y) \right)
```

and

```math
\tilde{\mathbf{v}} \left( r, \varphi \right)
    = \tilde{v}_r \left( r \right)
    = \mathbf{v} \left( \Phi^{-1}(x, y) \right)
```

the pressure and velocity distributions expressed in polar coordinates. The first identity
in the two above equations originates from the rotational symmetry of the problem and the
resulting independence of pressure and velocity on $`\varphi`$ in polar coordinates. Thus, in
polar coordinates we can write the mass balance equation as:

```math
\phi \frac{\partial \varrho}{\partial t}
   - \frac{\partial}{\partial r} \left( \frac{k}{\mu} \frac{\partial \tilde{p}}{\partial r} \right)
   = 0.
```

# Discretization

We employ a finite-volume scheme to spatially discretize the mass balance equation shown above.
The discrete equation describing mass conservation inside a control volume $`K`$ is obtained
by integration and reads:

```math
    | K | \left( \phi \, \partial \varrho / \partial t \right)_K
    + \sum_{\sigma \in \mathcal{S}_K} | \sigma | \left( \varrho v_r \right)_\sigma
    = 0,
```

where $`\sigma`$ are the faces of the control volume such that
$`\bigcup_{\sigma \in \mathcal{S}_K} \sigma \equiv \partial K`$ and where the notation $`( \cdot )_K`$
and $`( \cdot )_\sigma`$ was used to denote quantities evaluated for the control volume $`K`$ or a
face $`\sigma`$, respectively. The volume of the control volume is denoted with $`| K |`$ and
$`| \sigma |`$ is the area of a face.

Integration over polar coordinates requires taking into account the Jacobian determinant of the
coordinate transformation from polar to cartesian coordinates (see e.g.
[wikipedia.org](https://en.wikipedia.org/wiki/Polar_coordinate_system#Generalization)).
Let us discretize the domain by the intervals
$`K_i = (i\Delta r, (i+1)\Delta r) \times (0, 2 \Pi)`$,
$`i \in \{1, \dots, N \}`$, as control volumes.
As a result, their volumes are

```math
| K_i | = \Pi \left( ((i+1)*\Delta r)^2 - (i*\Delta r)^2 \right)
```

and the area of a face $`\sigma \in \mathcal{S}_{K_i}`$ is

```math
| \sigma | = 2 \Pi r_\sigma,
```

where $`r_\sigma`$ is the radius at which the face is defined.


### Header guard
This file contains the __problem class__ which defines the initial and boundary
conditions for the single-phase flow simulation.

### Include files
This header contains the porous medium problem class that this class is derived from:

```cpp
#include <dumux/porousmediumflow/problem.hh>
```

This header contains a convience function to calculate L2 errors

```cpp
#include <dumux/common/integrate.hh>
```

This header contains the class that specifies all spatially variable parameters
related to this problem.

```cpp
#include "spatialparams.hh"
```

### The problem class
We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
As this is a porous medium flow problem, we inherit from the base class `PorousMediumFlowProblem`.

```cpp
namespace Dumux {

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
```

We use convenient declarations that we derive from the property system.

```cpp
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
```

This is the constructor of our problem class:

```cpp
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), q_(0.0), pExact_(gridGeometry->numDofs())
    {
        k_ = getParam<Scalar>("SpatialParams.Permeability");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");
        q_ = getParam<Scalar>("Problem.Source");
        pW_ = getParam<Scalar>("Problem.WellPressure");
        rW_ = gridGeometry->bBoxMin()[0];

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                pExact_[scv.dofIndex()] = exactSolution(scv.dofPosition());
            }
        }
    }
```

First, we define the type of boundary conditions depending on the location. Two types of boundary  conditions
can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.

```cpp
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
```

We specify Dirichlet boundaries everywhere

```cpp
        values.setAllDirichlet();

        return values;
    }
```

Second, we specify the values for the Dirichlet boundaries. We need to fix values of our primary variable

```cpp
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
```

The exact solution values are set as Dirichlet values

```cpp
        return exactSolution(globalPos);
    }
```

We need to specify a constant temperature for our isothermal problem.
Fluid properties that depend on temperature will be calculated with this value.

```cpp
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }
```

This method add the exact pressure values to the vtk output

```cpp
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk) const
    {
        vtk.addField(pExact_, "pExact");
    }
```

The L2 error between the exact and numerical solution is calculated using this function,
using a specific order for the quadrature rule.

```cpp
    template<class SolutionVector>
    Scalar calculateL2Error(const SolutionVector& curSol, const int order)
    {
        return integrateL2Error(this->gridGeometry(), curSol, pExact_, order);
    }

private:
```

This function defines the exact pressure solution

```cpp
    PrimaryVariables exactSolution(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        const auto r = globalPos[0];
        priVars[0] = pW_ - 1.0/(2*M_PI)*nu_/k_*q_*std::log(r/rW_);
        return priVars;
    }

    Scalar q_, k_, nu_, rW_;
    GlobalPosition pW_;
    static constexpr Scalar eps_ = 1.5e-7;
    SolutionVector pExact_;
```

This is everything the one phase rotation symmetry problem class contains.

```cpp
};
```

We leave the namespace Dumux.

```cpp
} // end namespace Dumux
```





```cpp



#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar,
                             OnePTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar,
                                           OnePTestSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;
    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

private:
    Scalar permeability_;
};

} // end namespace Dumux
```




We look now at the main file for the tracer problem. We set up two problems in this file and solve them sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. We use it to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes, we calculate the transport of a tracer in the following tracer problem.
### Includes

```cpp
#include <config.h>
```

This includes the `TypeTags` and properties to be used for the single-phase rotation symmetry example.

```cpp
#include "properties.hh"
```

Further, we include a standard header file for C++, to get time and date information

```cpp
#include <ctime>
```

and another one for in- and output.

```cpp
#include <iostream>
```

Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers.
Here, we include classes related to parallel computations, time measurements and file I/O.

```cpp
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
```

In Dumux, the property system is used to specify classes and compile-time options to be used by the model.
For this, different properties are defined containing type definitions, values and methods.
All properties are declared in the file `properties.hh`.

```cpp
#include <dumux/common/properties.hh>
```

The following file contains the parameter class, which manages the definition and retrieval of input
parameters by a default value, the inputfile or the command line.

```cpp
#include <dumux/common/parameters.hh>
```

The file `dumuxmessage.hh` contains the class defining the start and end message of the simulation.

```cpp
#include <dumux/common/dumuxmessage.hh>
```

The following file contains the class, which defines the sequential linear solver backends.

```cpp
#include <dumux/linear/seqsolverbackend.hh>
```

The following linear pde solver allows to solve linear equations

```cpp
#include <dumux/linear/pdesolver.hh>
```

Further we include the assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation).

```cpp
#include <dumux/assembly/fvassembler.hh>
```

The containing class in the following file defines the different differentiation methods used to compute the derivatives of the residual.

```cpp
#include <dumux/assembly/diffmethod.hh>
```

We need the following class to simplify the writing of dumux simulation data to VTK format.

```cpp
#include <dumux/io/vtkoutputmodule.hh>
```

The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.

```cpp
#include <dumux/io/grid/gridmanager.hh>
```

### Beginning of the main function

```cpp
int main(int argc, char** argv) try
{
    using namespace Dumux;
```

Convenience aliases for the type tag of the problem.

```cpp
    using TypeTag = Properties::TTag::OnePRotSymBox;
```

We initialize MPI. Finalization is done automatically on exit.

```cpp
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
```

We print the dumux start message.

```cpp
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);
```

We parse the command line arguments.

```cpp
    Parameters::init(argc, argv);
```

### Create the grid
The `GridManager` class creates the grid from information given in the input file.
This can either be a grid file, or in the case of structured grids, by specifying the coordinates
of the corners of the grid and the number of cells to be used to discretize each spatial direction.
Here, we solve both the single-phase and the tracer problem on the same grid.
Hence, the grid is only created once using the grid type defined by the type tag of the 1p problem.

```cpp
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
```

start timer

```cpp
    Dune::Timer timer;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    int numRefinements = getParam<int>("Grid.RefinementSteps");
    constexpr bool isBox = GridGeometry::discMethod == Dumux::DiscretizationMethod::box;
    const int orderQuadratureRule = isBox ? 3 : 1;
    std::vector<Scalar> l2Errors(numRefinements, 0.0);
    for(int i = 0; i < numRefinements; i++)
    {
```

We compute on the leaf grid view.

```cpp
        const auto& leafGridView = gridManager.grid().leafGridView();
```

### Set-up and solving of the 1p rotation symmetry problem
In the following section, we set up and solve the problem. As the result of this problem, we obtain the pressure distribution in the domain.
#### Set-up
We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.

```cpp
        auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
        gridGeometry->update();
```

In the problem, we define the boundary and initial conditions.

```cpp
        using Problem = GetPropType<TypeTag, Properties::Problem>;
        auto problem = std::make_shared<Problem>(gridGeometry);
```

The solution vector containing the pressure values

```cpp
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        SolutionVector x(gridGeometry->numDofs());
```

The grid variables store variables (primary and secondary variables) on sub-control volumes and faces (volume and flux variables).

```cpp
        using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);
```

#### Assembling the linear system
We create the assembler.

```cpp
        using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
```

#### Solution
We set the linear solver "UMFPack" as the linear solver. Afterwards we solve the linear system.
The linear system of equations is solved

```cpp
        using LinearSolver = UMFPackBackend;
        auto linearSolver = std::make_shared<LinearSolver>();
        LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
        solver.setVerbose(0);
        solver.solve(x);
```

#### Post-processing and output
We calculate the L2 errors using the numerical solution

```cpp
        l2Errors[i] = problem->calculateL2Error(x, orderQuadratureRule);
        std::cout.precision(8);
        std::cout << "L2 error for "
                  << std::setw(6) << gridGeometry->numDofs()
                  << " dofs: "
                  << std::scientific
                  << l2Errors[i]
                  << " rate: "
                  << ((i>0) ? std::to_string(std::log(l2Errors[i]/l2Errors[i-1])/(std::log(0.5))) : " - ")
                  << std::endl;
```

The vtk file is written on the finest grid.

```cpp
        if(i == numRefinements-1)
        {
```

output result to vtk

```cpp
            VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
            using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
            vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
            using IOFields = GetPropType<TypeTag, Properties::IOFields>;
            IOFields::initOutputModule(vtkWriter); // Add model specific output fields
            problem->addVtkFields(vtkWriter); //!< Add problem specific output fields
            vtkWriter.write(0.0);
        }
```

Globally refine the grid and repeate the solution procedure

```cpp
        gridManager.grid().globalRefine(1);
    }
```

We stop the timer and display the total time of the simulation as well as the cumulative CPU time.

```cpp
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";
```

### Final Output, print parameters

```cpp
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
}
```

### Exception handling
In this part of the main file we catch and print possible exceptions that could
occur during the simulation.

```cpp
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
```


# TODO: Paraview stuff