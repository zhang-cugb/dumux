// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Assembly
 * \brief Assembler class for residuals and jacobian matrices for
 *        grid-based numerical schemes and stationary problems.
 */
#ifndef DUMUX_STATIONARY_ASSEMBLER_HH
#define DUMUX_STATIONARY_ASSEMBLER_HH

#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/assembly/partialreassembler.hh>

// TODO: include box/cc local assemblers
#include "fem/stationarylocalassembler.hh"
#include "jacobianpattern.hh"

namespace Dumux {
namespace Impl {

    template<class Assembler, DiscretizationMethod dm> struct LocalAssemblerChooser;

    // TODO: Currently the fv local residuals require type tag etc...
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::tpfa> { using type = CCLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::mpfa> { using type = CCLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::box> { using type = BoxLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::staggered> { using type = ... };
    template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::fem> { using type = FEStationaryLocalAssembler<Assembler>; };

    template<class Assembler, DiscretizationMethod dm>
    using LocalAssemblerType = typename LocalAssemblerChooser<Assembler, dm>::type;

} // end namespace detail

//! Default types used for the linear system
template<class Scalar, int numEq>
struct DefaultLinearSystemTraits
{
private:
    using PrimaryVariables = Dune::FieldVector<Scalar, numEq>;
    using BlockType = Dune::FieldMatrix<Scalar, numEq, numEq>;

public:
    using ResidualVector = Dune::BlockVector<PrimaryVariables>;
    using JacobianMatrix = Dune::BCRSMatrix<BlockType>;
};

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for grid-based numerical schemes
 * \tparam LO The local operator (evaluation of the terms of the equation)
 * \tparam diffMethod The differentiation method to compute derivatives
 * \tparam LST The linear system traits (types used for jacobians and residuals)
 */
template< class LO, DiffMethod diffMethod,
          class LST = DefaultLinearSystemTraits<typename LO::GridVariables::Scalar,
                                                LO::GridVariables::PrimaryVariables::size()> >
class StationaryAssembler
{
    using ThisType = StationaryAssembler<LO, diffMethod, LST>;
    using GridView = typename LO::GridVariables::GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the types used for the linear system
    using Scalar = typename LO::GridVariables::Scalar;
    using JacobianMatrix = typename LST::JacobianMatrix;
    using ResidualVector = typename LST::ResidualVector;
    using ResidualType = ResidualVector; // Required by NewtonSolver etc. (TODO: Which name is better!?)

    //! export the local operator type
    using LocalOperator = LO;

    //! the local operator states the type of variables needed for evaluation
    using GridVariables = typename LO::GridVariables;

    //! export a grid-independent alias such that e.g. the newton solver works
    //! for grid-based schemes as well as algebraic non-linear systems of equations
    using Variables = GridVariables;

    //! The variables depend on a user-defined problem (boundary conditions/parameters)
    using Problem = typename GridVariables::Problem;

    //! This assembler is for grid-based schemes, so the variables live on a grid geometry
    using GridGeometry = typename GridVariables::GridGeometry;

    //! The local assembler is the local view on this assembler class
    using LocalView = Impl::LocalAssemblerType<ThisType, GridGeometry::discMethod>;

    //! Export type used for solution vectors (the variables depend upon the solution)
    using SolutionVector = typename GridVariables::SolutionVector;

    /*!
     * \brief The Constructor from a grid geometry.
     * \param gridGeometry A grid geometry instance
     * \note This assembler class is, after construction, defined for a specific equation
     *       (given by the template argument of the LocalOperator) and a specific grid
     *       geometry - which defines the connectivity of the degrees of freedom of the
     *       underlying discretization scheme on a particular grid. The evaluation point,
     *       consisting of a particular solution/variables/parameters may vary, an therefore,
     *       an instance of the grid variables is passed to the assembly functions.
     */
    StationaryAssembler(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    {}

    /*!
     * \brief Assembles the Jacobian matrix and the residual around the given evaluation point
     *        which is determined by the grid variables, containing all quantities required
     *        to evaluate the equations to be assembled.
     * \param gridVariables The variables corresponding to the given solution state
     * \note We assume the grid geometry on which the grid variables are defined
     *       to be the same as the one used to instantiate this class
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        resetJacobian_(partialReassembler);
        resetResidual_();

        assemble_([&](const Element& element)
        {
            auto localAssembler = localView(*this);
            localAssembler.bind(element, gridVariables);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, partialReassembler);
        });

        enforcePeriodicConstraints_(*jacobian_, *residual_);
    }

    /*!
     * \brief Assembles the Jacobian matrix of the discrete system of equations
     *        around a given state represented by the grid variables object.
     */
    void assembleJacobian(const GridVariables& gridVariables)
    {
        resetJacobian_();
        assemble_([&](const Element& element)
        {
            auto localAssembler = localView(*this);
            localAssembler.bind(element, gridVariables);
            localAssembler.assembleJacobianAndResidual(*jacobian_);
        });
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the internal residual vector to store the result.
     */
    void assembleResidual(const GridVariables& gridVariables)
    {
        resetResidual_();
        assembleResidual(*residual_, gridVariables);
    }

    /*!
     * \brief Assembles the residual for a given state represented by the provided
     *        grid variables object, using the provided residual vector to store the result.
     */
    void assembleResidual(ResidualVector& r, const GridVariables& gridVariables) const
    {
        // TODO: I leave the line below from the old assembler (commented) because
        //       I think it illustrates how the solution and the grid variables are
        //       coupled and that it is probably better to make the assembly routines
        //       functions of grid variables instead of solutions.
        //
        // update the grid variables for the case of active caching
        // gridVariables_->update(curSol);

        assemble_([&](const Element& element)
        {
            auto localAssembler = localView(*this);
            localAssembler.bind(element, gridVariables);
            localAssembler.assembleResidual(r);
        });
    }

    //! TODO: This is not so nice! But, it allows the newton to be completely agnostic
    //!       about how to update the variables after a new iterate has been computed
    //!       Until now, this was the updateGridVariables() interface in FVAssembler
    //!       Is there a nicer way to generalize the update such that the newton can
    //!       stay agnostic about it but we don't need to have it in the assembler!?
    //!       One thing to consider is that for time-dependent problems, the grid variables
    //!       are possibly time-dependent and need more than just the current solution.
    //!       The InstationaryAssembler would probably carry the time information and could
    //!       do the update... But still, maybe there is a better place...
    //!
    //! \brief Update a grid variables instance to the given solution
    void update(GridVariables& gridVariables, const SolutionVector& x)
    {
        gridVariables.update(x);
    }

    //! TODO: Do we want to remove this interface?
    //!       Should it be the assembler's job to compute the norm?
    //! compute the residual and return it's vector norm
    Scalar residualNorm(const GridVariables& gridVars) const
    {
        ResidualVector residual(numDofs());
        assembleResidual(residual, gridVars);

        // for box communicate the residual with the neighboring processes
        // TODO: INSTEAD OF ISBOX, DETERMINE GENERALLY IF THE SCHEME NEEDS THIS AND UNCOMMENT
        // if (isBox && gridView().comm().size() > 1)
        // {
        //     using VertexMapper = typename GridGeometry::VertexMapper;
        //     VectorCommDataHandleSum<VertexMapper, SolutionVector, GridGeometry::GridView::dimension>
        //         sumResidualHandle(gridGeometry_->vertexMapper(), residual);
        //     gridView().communicate(sumResidualHandle,
        //                            Dune::InteriorBorder_InteriorBorder_Interface,
        //                            Dune::ForwardCommunication);
        // }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        // check and/or set the BCRS matrix's build mode
        if (jacobian_->buildMode() == JacobianMatrix::BuildMode::unknown)
            jacobian_->setBuildMode(JacobianMatrix::random);
        else if (jacobian_->buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        jacobian_->setBuildMode(JacobianMatrix::random);
        residual_ = std::make_shared<SolutionVector>();

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // create occupation pattern of the jacobian
        // TODO: HOW TO DETERMINE PATTERN DEPENDING ON TIME SCHEME??
        const auto occupationPattern = getJacobianPattern</*isImplicit*/true>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return gridGeometry().numDofs(); }

    //! The global finite volume geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The gridview
    const GridView& gridView() const
    { return gridGeometry().gridView(); }

    //! The jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! The residual vector (rhs)
    ResidualVector& residual()
    { return *residual_; }

protected:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if (!residual_)
        {
            residual_ = std::make_shared<ResidualVector>();
            setResidualSize();
        }

        (*residual_) = 0.0;
    }

    // reset the Jacobian matrix to 0.0
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if (!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        if (partialReassembler)
            partialReassembler->resetJacobian(*this);
        else
            *jacobian_ = 0.0;
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     */
    template<typename AssembleElementFunc>
    void assemble_(AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                assembleElement(element);

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView().comm().size() > 1)
            succeeded = gridView().comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    void enforcePeriodicConstraints_(JacobianMatrix& jac, SolutionVector& res)
    { /*TODO: Implement*/ }

private:
    //! the grid geometry on which it is assembled
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualVector> residual_;
};

} // namespace Dumux

#endif
