#include <config.h>

#include <iostream>
#include <cmath>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dumux/nonlinear/newtonsolver.hh>

/*

  This test currently solves a scalar non-linear equation using the
  Dumux::NewtonSolver. The Mock classes expose which dependencies the
  current implementation has on different other classes it interacts with.
  Several dependencies seem unnecessary. In particular, the current
  implementation is basically hard-coded to double indexable residual vectors
  with several math operators. A good idea would seem to somehow delegate
  this dependency to something like a linear algebra backend or at least
  the assembler. The assembler requires a lot of interface functions which
  are not always needed. The linear solver interdependency is much better (small).

  This test is to ensure that the dependencies do not grow more in the future.

 */

namespace Dumux {

class MockScalarAssembler
{
public:
    using Variables = double;
    using ResidualType = double;
    using Scalar = double;
    using JacobianMatrix = double;

    void setLinearSystem() {}

    bool isStationaryProblem() { return false; }

    ResidualType prevSol() { return ResidualType(0.0); }

    void resetTimeStep(const ResidualType& x) {}

    void assembleResidual(const ResidualType& x)
    {
        res_ = x*x - 5.0;
    }

    void assembleResidual(ResidualType& r, const ResidualType& x) const
    {
        r = x*x - 5.0;
    }

    void assembleJacobianAndResidual (const ResidualType& x)
    {
        assembleResidual(x);
        jac_ = 2.0*x;
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

    void updateGridVariables(const ResidualType& sol) {}

    Variables& gridVariables() { return vars_; }

private:
    JacobianMatrix jac_;
    ResidualType res_;
    Variables vars_;
};

class MockScalarLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    bool solve (const double& A, double& x, const double& b) const
    { x = b/A; return true; }

    double norm (const double& r)
    { return std::abs(r); }
};

} // end namespace Dumux

int main(int argc, char* argv[]) try
{
    using namespace Dumux;

    // maybe initialize MPI
    Dune::MPIHelper::instance(argc, argv);

    // initialize  parameters
    // TODO this is necessary because there are some global default used in the Newton solver
    // Do we really need them to be global defaults???
    Parameters::init(argc, argv);

    // use the Newton solver to find a solution to a scalar equation
    using Assembler = MockScalarAssembler;
    using LinearSolver = MockScalarLinearSolver;
    using Solver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>();
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    double initialGuess = 0.1;
    double x(initialGuess);

    std::cout << "Solving: x^2 - 5 = 0" << std::endl;
    solver->solve(x);
    std::cout << "Solution: " << std::setprecision(15) << x
              << ", exact: " << std::sqrt(5.0)
              << ", error: " << std::abs(x-std::sqrt(5.0))/std::sqrt(5.0)*100 << "%" << std::endl;

    if (Dune::FloatCmp::ne(x, std::sqrt(5.0), 1e-13))
        DUNE_THROW(Dune::Exception, "Didn't find correct root: " << std::setprecision(15) << x << ", exact: " << std::sqrt(5.0));

    return 0;

}
catch (const Dune::Exception& e)
{
    std::cout << e << std::endl;
    return 1;
}
catch (...)
{
    std::cout << "Unknown exception thrown!" << std::endl;
    return 1;
}
