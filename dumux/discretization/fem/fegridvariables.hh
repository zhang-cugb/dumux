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
 * \ingroup Discretization
 * \brief The grid variable class for finite element schemes.
 * \note This default implementation does not store any additional data on the grid.
 */
#ifndef DUMUX_FE_GRID_VARIABLES_HH
#define DUMUX_FE_GRID_VARIABLES_HH

#include <cassert>
#include <dumux/discretization/localview.hh>
#include <dumux/timestepping/timelevel.hh>

#include "elementsolution.hh"
#include "ipvariablesbase.hh"

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Local view on the grid variables for finite element schemes.
 * \note This default implementation only stores a pointer on the grid variables.
 * \tparam GV The grid variables class
 */
template<class GV>
class FEGridVariablesLocalView
{
    using GridGeometry = typename GV::GridGeometry;
    using FEElementGeometry = typename GridGeometry::LocalView;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using PrimaryVariables = typename GV::PrimaryVariables;
    using ElementSolution = FEElementSolution<FEElementGeometry, PrimaryVariables>;

public:
    using GridVariables = GV;

    //! The constructor
    FEGridVariablesLocalView(const GridVariables& gridVariables)
    : gridVariables_(&gridVariables)
    {}

    /*!
     * \brief Bind this local view to a grid element.
     * \param element The grid element
     * \param feGeometry Local view on the grid geometry
     * \todo TODO: In fv schemes here the current solution enters
     *             (in the vol vars/flux vars cache) in order to bind
     *             the local views to a potentially deflected solution
     *             Do we require the grid variables to carry the notion
     *             of the deflected solution or should this be injected here!?
     */
    void bind(const Element& element,
              const FEElementGeometry& feGeometry)
    {
        const auto& x = gridVariables().dofs();
        const auto& gg = gridVariables().gridGeometry();
        elemSol_ = elementSolution(element, x, gg);
    }

    /*!
     * \brief Return the element solution, i.e. the element-local
     *        view on the solution vector.
     */
    const ElementSolution& elemSol() const
    { return elemSol_; }

    /*!
     * \todo TODO: This is currently required such that the local assembler
     *             can deflect the solution in here in order to compute numerical
     *             derivatives. This means that the state of this local object does
     *             not correspond to the global object!?
     */
    ElementSolution& elemSol()
    { return elemSol_; }

    /*!
     * \brief Return a reference to the grid variables.
     */
    const GridVariables& gridVariables() const
    { return *gridVariables_; }

private:
    const GridVariables* gridVariables_;
    ElementSolution elemSol_;
};

//! TODO: have these defined at a central place
template<class Problem>
struct ProblemTraits;

/*!
 * \ingroup Discretization
 * \brief The grid variable class for finite element schemes.
 * \note This default implementation does not store any additional data on the grid.
 * \tparam P the problem to be solved
 * \tparam X the type used to represent solution vectors on the grid
 * \tparam IPV the variable type to be constructed at integration points
 */
template<class P, class X, class IPV = IntegrationPointVariablesBase<typename X::value_type>>
class FEGridVariables
{
    using ThisType = FEGridVariables<P, X, IPV>;

public:
    //! export the underlying problem
    using Problem = P;

    //! export the type used for solution vectors
    using SolutionVector = X;

    //! export type of the finite volume grid geometry
    using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;

    //! export primary variable type
    using PrimaryVariables = typename SolutionVector::value_type;

    //! export scalar type
    using Scalar = typename PrimaryVariables::value_type;

    //! export type of variables constructed at integration points
    using IntegrationPointVariables = IPV;

    //! The local view on these grid variables
    using LocalView = FEGridVariablesLocalView<ThisType>;

    //! Type representing a time level
    using TimeLevel = Dumux::TimeLevel<Scalar>;

    //! constructor
    FEGridVariables(std::shared_ptr<const Problem> problem,
                    const TimeLevel& timeLevel = TimeLevel{0.0})
    : problem_(problem)
    , timeLevel_(timeLevel)
    {
        problem->applyInitialSolution(x_);
    }

    //! TODO
    template<class ApplyInitialSol>
    FEGridVariables(std::shared_ptr<const Problem> problem,
                    const ApplyInitialSol& applyInitialSol,
                    const TimeLevel& timeLevel = TimeLevel{0.0})
    : problem_(problem)
    , timeLevel_(timeLevel)
    {
        applyInitialSol(x_);
    }

    //! return the underlying problem
    const Problem& problem() const
    { return *problem_; }

    //! return the finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return problem_->gridGeometry(); }

    //! return the solution for which the grid variables were updated
    const SolutionVector& dofs() const
    { return x_; }

    //! return the time level the grid variables represent
    const TimeLevel& timeLevel() const
    { return timeLevel_; }

    //! update all variables subject to a new solution
    void updateDofs(const SolutionVector& x)
    {
        x_ = x;
    }

    //! TODO Doc me
    void updateTime(const TimeLevel& timeLevel)
    {
        timeLevel_ = timeLevel;
    }

    //! update all variables subject to a new solution
    void update(const SolutionVector& x, const TimeLevel& timeLevel)
    {
        updateDofs(x);
        updateTime(timeLevel);
    }

private:
    std::shared_ptr<const Problem> problem_; //!< pointer to the problem to be solved
    SolutionVector x_;                       //!< the solution corresponding to this class' current state
    TimeLevel timeLevel_;                    //!< info during time integration for instationary problems
};

} // end namespace Dumux

#endif
