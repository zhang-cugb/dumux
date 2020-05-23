#include <config.h>

#include <iostream>
#include <iomanip>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/numericdifferentiation.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/timestepping/multistagetimestepper.hh>
#include <dumux/timestepping/multistagelocaloperator.hh>
#include <dumux/timestepping/multistagemethods.hh>

/*
  This test currently solves a scalar non-linear and time-dependent
  equation using the Dumux::NewtonSolver and Dumux::MultiStageTimeStepper.
  TODO: DOC THIS...
 */

namespace Dumux {

// TEMPORARY COMMENT FOR UNDERSTANDING:
// In the general case, we expect local residuals (the actual equation) to
// depend on variables, whose state depend on a current solution (primary variables)
// and potentially the time as well as more user input (e.g. via a problem). They are
// required to have an update function that allows updating all variables subject to a
// new solution. For instationary problems, the variables also depend on time. Here, we
// define the interface to receive the time, but it could also receive the stage parameters!?
//
// This class could/should be defined more slim for this test. The data types could be
// simply defined globally for this test. However, to mimick a situation more similar to
// that of grid based schemes for now to identify the couplings.
template<class SV>
class MockVariables
{
public:
    //! export solution vector type
    using SolutionVector = SV;
    //! export primary variables type (NOT USED IN THIS EXAMPLE)
    using PrimaryVariables = typename SV::value_type;
    //! export scalar type (ALSO NOT USED - ONLY WITHIN THIS CLASS)
    using Scalar = typename PrimaryVariables::value_type;

    //! constructor
    MockVariables(const SolutionVector& curSol,
                  Scalar time)
    : curSol_(&curSol)
    , time_(time)
    {}

    //! update variables to a new solution
    void update(const SolutionVector& sol, Scalar time)
    {
        curSol_ = &sol;
        time_ = time;
    }

    //! current solution the variables are updated to
    const SolutionVector& curSol() const
    { return *curSol_; }

    //! time level of these variables
    Scalar time() const
    { return time_; }

private:
    const SolutionVector* curSol_;
    Scalar time_;
};

// TEMPORARY COMMENT FOR UNDERSTANDING:
// A local residual alows for (local) evaluation of the equation to be solved.
// We assume (for instationary problems) equations of the form:
// ds(u, t)/dt + F(u, t) = 0,
// where F(u, t) can be composed of fluxes and sources, but no time-derivatives w.r.t
// the unknown u are contained. These are collected in the term ds(u, t)/dt.
// The local residuals must have evalSpatial() and evalTemporal() to evaluate the terms
// that DO NOT appear or DO appear in time derivatives. Thus, evalTemporal() evaluates
// s(u, t), while evalSpatial() evaluates F(u, t).
// This residual implementation evaluates du/dt - e^t = 0.
template<class Vars>
class MockLocalOperator
{
public:
    //! export the type containing the variables required to evaluate the residual
    using Variables = Vars;
    //! export the underlying scalar type
    using Scalar = typename Variables::Scalar;
    //! export the type representing an evaluated local residual
    using Residual = Dune::FieldVector<Scalar, 1>;

    //! Construct a local residual
    //! In general, local residuals are instantiated with local views
    //! on variable objects. In this simple setting, our variables consist
    //! of two single scalar which is the solution and the time. Thus, we simply
    //! construct it from that. Construction is done in our test assembler class.
    //! In grid based methods, there is a local assembler between the assembler
    //! and the local residual, which could construct it from a local view.
    MockLocalOperator(const Variables& variables)
    : variables_(variables)
    {}

    // evaluate the term appearing in time derivative
    Residual evalStorage() const
    {
        // our time-dependent term is simply u, stored in the solution vector
        const auto u = variables_.curSol()[0][0];
        return u;
    }

    // evaluate the term NOT appearing in time derivative
    Residual evalFluxesAndSources() const
    {
        const auto t = variables_.time();
        return -1.0*std::exp(t);
    }

    //! Return an empty residual
    //! TODO: This is a newly introduced interface. Residuals (for grid-based
    //!       schemes ElementResidualVector) are generally resizable vectors,
    //!       containing the residual at all dofs of an element. In instationary
    //!       settings using MultiStageLocalOperator, the MultiStageLocalOperator
    //!       computes the residual by summing up several weighted terms. It is
    //!       therefore handy to be able to get a (correctly resized if needed)
    //!       empty residual in which to add up things. Alternatively, we could
    //!       use a resizable vector type here as return type. In case we want
    //!       such a functionality, we could also put this somewhere else. Maybe
    //!       a free function like to call like
    //!       auto emptyRes = getEmptyResidual(localResidual);
    Residual getEmptyResidual() const
    { return 0.0; }

private:
    const Variables& variables_;
};

template<class LocalOperator>
class InstationaryMockAssembler
{
public:
    //! export the variable type used for residual evaluations
    using Variables = typename LocalOperator::Variables;

    //! export underlying data types
    using Scalar = typename Variables::Scalar;
    using SolutionVector = typename Variables::SolutionVector;

    //! export linear sytem types
    using JacobianMatrix = Scalar;
    using ResidualType = SolutionVector;

    //! export parameters required for time integration
    using StageParams = MultiStageParams<Scalar>;

    ///////////////////////////////////////////////
    // existing assembler interfaces
    void setLinearSystem() {}

    void assembleResidual(const Variables& vars)
    {
        // wrap multi-stage residual around local operators
        std::vector<LocalOperator> localOperators;
        localOperators.reserve(stageParams_->size());
        for (std::size_t i = 0; i < stageParams_->size()-1; ++i)
            localOperators.emplace_back(stageVariables_[i]);
        localOperators.emplace_back(vars);

        using MSLocalOperator = MultiStageLocalOperator<LocalOperator>;
        res_[0][0] = MSLocalOperator(localOperators, *stageParams_).evalLocalResidual();
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);

        // For this problem, the derivative is always 1.0
        jac_ = 1.0;

        // numerical differentiation also works (although
        // the code is a bit cumbersome)...
        //
        // // make a copy of the variables to allow deflection
        // auto nonConstVars = vars;
        // const auto origRes = res_;
        //
        // // lambda to evaluate the local residual
        // // as function of a primary variable
        // auto evalRes = [&] (auto x)
        // {
        //     // make solution vector from given solution
        //     // and update variables
        //     SolutionVector curX;
        //     curX[0][0] = x;
        //     nonConstVars.update(curX, vars.time());
        //
        //     // assemble and return residual
        //     this->assembleResidual(nonConstVars);
        //     return this->res_[0][0];
        // };
        //
        // // write the partial derivative in jac_
        // NumericDifferentiation::partialDerivative(evalRes, vars.curSol()[0][0], jac_, origRes[0][0]);
        //
        // // restore original state of the residual
        // res_ = origRes;
    }

    JacobianMatrix& jacobian() { return jac_; }
    ResidualType& residual() { return res_; }
    ///////////////////////////////////////////////

    ///////////////////////////////////////////////
    // THE FOLLOWING ARE ADDITIONAL NEW INTERFACES
    // - MAYBE WE CAN USE THESE FOR INSTATIONARY ASSEMBLERS

    // prepare for a stage in a time integration step
    void prepareStage(SolutionVector& x,
                      Variables& vars,
                      std::shared_ptr<const StageParams> params)
    {
        stageParams_ = params;
        const auto curStage = params->size() - 1;

        // clean up if a new time step is started
        if (curStage == 0)
            stageVariables_.clear();

        // make sure the variables live on the correct time level
        // (the solution in the vars should be the right one already, but so far
        // we have no possibility to just update the time)
        vars.update(x, params->timeAtStage(curStage));
    }

    // register the state of a finished stage
    void registerStage(SolutionVector& x,
                       Variables& vars)
    {
        stageVariables_.emplace_back(vars);
    }

    //! update variables to new solution
    void update(Variables& vars, const SolutionVector& x)
    {
        const auto curStage = stageParams_->size() - 1;
        vars.update(x, stageParams_->timeAtStage(curStage));
    }

    ///////////////////////////////////////////
    //! TODO: INTERFACES THAT WE WANT TO REMOVE?
    bool isStationaryProblem() { return false; }
    ResidualType prevSol() { return ResidualType(0.0); }
    void resetTimeStep(const ResidualType& sol) {}
    void updateGridVariables(const SolutionVector& sol) {}
    double residualNorm(const Variables& vars)
    {
        assembleResidual(vars);
        return res_[0][0];
    }
    ////////////////////////////////////////

private:
    JacobianMatrix jac_;
    ResidualType res_;

    std::shared_ptr<const StageParams> stageParams_;
    std::vector<Variables> stageVariables_;
};

template<class Scalar>
class MockLinearSolver
{
public:
    void setResidualReduction(Scalar residualReduction) {}

    template<class Vector>
    bool solve (const Scalar& A, Vector& x, const Vector& b) const
    {
        x[0][0] = b[0][0]/A;
        return true;
    }
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

    using Scalar = double;
    using PrimaryVariables = Dune::FieldVector<Scalar, 1>;
    using SolutionVector = Dune::FieldVector<PrimaryVariables, 1>;

    using Variables = MockVariables<SolutionVector>;
    using LocalOperator = MockLocalOperator<Variables>;

    using Assembler = InstationaryMockAssembler<LocalOperator>;
    using LinearSolver = MockLinearSolver<Scalar>;
    using Solver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

    // solution vector
    using SolutionVector = typename Assembler::SolutionVector;
    SolutionVector x = 0.0;

    // variables storing everything required to evaluate the equations
    using Variables = typename Assembler::Variables;
    Variables vars(x, /*time=*/0.0);

    // create PDE solver
    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>();
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    Scalar dt = getParam<Scalar>("TimeLoop.Dt");

    // 1. Do explicit Euler time integration
    auto expEuler = std::make_shared< MultiStage::ExplicitEuler<Scalar> >();
    MultiStageTimeStepper<Solver> timeStepper(solver, expEuler);

    std::cout << "Solving 5 time steps with " << expEuler->name() << std::endl;
    std::array<Scalar, 5> expEulerValues;
    for (unsigned int i = 0; i < 5; ++i)
    {
        timeStepper.step(x, vars, i*dt, dt);
        expEulerValues[i] = x[0][0];
    }

    // 2. Do theta scheme with 0.5 time integration
    x = 0.0;
    auto theta = std::make_shared< MultiStage::Theta<Scalar> >(0.5);
    timeStepper.setMethod(theta);

    std::cout << "Solving 5 time steps with " << theta->name() << std::endl;
    std::array<Scalar, 5> thetaValues;
    for (unsigned int i = 0; i < 5; ++i)
    {
        timeStepper.step(x, vars, i*dt, dt);
        thetaValues[i] = x[0][0];
    }

    // 3. Do implicit Euler time integration
    x = 0.0;
    auto impEuler = std::make_shared< MultiStage::ImplicitEuler<Scalar> >();
    timeStepper.setMethod(impEuler);

    std::cout << "Solving 5 time steps with " << impEuler->name() << std::endl;
    std::array<Scalar, 5> impEulerValues;
    for (unsigned int i = 0; i < 5; ++i)
    {
        timeStepper.step(x, vars, i*dt, dt);
        impEulerValues[i] = x[0][0];
    }

    // 4. Do explicit Runge-Kutta 4th order time integration
    x = 0.0;
    auto expRK4 = std::make_shared< MultiStage::RungeKuttaExplicitFourthOrder<Scalar> >();
    timeStepper.setMethod(expRK4);

    std::cout << "Solving 5 time steps with " << expRK4->name() << std::endl;
    std::array<Scalar, 5> expRK4Values;
    for (unsigned int i = 0; i < 5; ++i)
    {
        timeStepper.step(x, vars, i*dt, dt);
        expRK4Values[i] = x[0][0];
    }

    // print computed values vs exact solution
    auto exact = [] (auto t) { return std::exp(t) - 1; };
    std::cout << "\n\n"
              << "Exact solutions:              ";
    for (unsigned int i = 1; i < 6; i++)
        std::cout << std::left << std::setw(9) << std::setfill(' ') << exact(i*dt) << (i < 5 ? ", " : "\n");

    auto printSchemeSol = [] (const auto& schemeValues)
    {
        for (unsigned int i = 0; i < schemeValues.size(); ++i)
            std::cout << std::left << std::setw(9) << std::setfill(' ') << schemeValues[i] << (i < schemeValues.size()-1 ? ", " : "\n");
    };

    std::cout << "Solutions for explicit Euler: "; printSchemeSol(expEulerValues);
    std::cout << "Solutions for theta scheme:   "; printSchemeSol(thetaValues);
    std::cout << "Solutions for implicit Euler: "; printSchemeSol(impEulerValues);
    std::cout << "Solutions for explicit RK4:   "; printSchemeSol(expRK4Values);

    // compute errors of the schemes via trapezoidal rule
    auto computeError = [&] (const auto& schemeValues)
    {
        Scalar error = 0.0;
        for (unsigned int i = 0; i < schemeValues.size(); ++i)
        {
            const auto exactSol = exact((i+1)*dt);
            error += (exactSol - schemeValues[i])*(exactSol - schemeValues[i])*dt;
        }

        return error;
    };

    const auto errorExpEuler = computeError(expEulerValues);
    const auto errorTheta = computeError(thetaValues);
    const auto errorImpEuler = computeError(impEulerValues);
    const auto errorExpRK4 = computeError(expRK4Values);
    std::cout << "\n\nErrors: \n"
              << "Explicit Euler:  " << errorExpEuler << "\n"
              << "Theta (0.5):     " << errorTheta << "\n"
              << "Implicit Euler:  " << errorImpEuler << "\n"
              << "Explicit RK4:    " << errorExpRK4 << "\n";

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
