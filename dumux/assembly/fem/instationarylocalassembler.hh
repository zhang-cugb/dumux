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
 * \brief An assembler for Jacobian and residual contribution
 *        per element (finite element method) in the context of
 *        instationary problems.
 */
#ifndef DUMUX_FE_INSTATIONARY_LOCAL_ASSEMBLER_HH
#define DUMUX_FE_INSTATIONARY_LOCAL_ASSEMBLER_HH

#include <cassert>
#include <memory>

#include <dune/grid/common/gridenums.hh>

#include <dumux/discretization/fem/elementsolution.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/timestepping/multistagetimestepper.hh>
#include <dumux/timestepping/multistagelocaloperator.hh>


namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution
 *        per element (finite element method) in the context of
 *        instationary problems.
 * \tparam Assembler The grid-wide assembler type
 * \todo TODO: This assumes Standard-Galerkin discretizations and non-composite function spaces.
 *             Could this be generalized to other spaces and Petrov-Galerkin?
 */
template<class Assembler, DiffMethod diffMethod>
class FEInstationaryLocalAssembler
{
    using SolutionVector = typename Assembler::SolutionVector;
    using GridVariables = typename Assembler::GridVariables;
    using GridGeometry = typename GridVariables::GridGeometry;
    using Problem = typename GridVariables::Problem;

    using FEElementGeometry = typename GridGeometry::LocalView;
    using ElementVariables = typename GridVariables::LocalView;

    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using ElementSolution = FEElementSolution<FEElementGeometry, PrimaryVariables>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numEq = PrimaryVariables::size();

    using StageParams = typename Assembler::StageParams;
    using LocalOperator = typename Assembler::LocalOperator;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ResidualVector = typename Assembler::ResidualVector;

    static_assert(diffMethod == DiffMethod::numeric, "Analytical assembly not implemented");

public:
    using ElementResidualVector = typename LocalOperator::ElementResidualVector;

    /*!
     * \brief Constructor.
     */
    explicit FEInstationaryLocalAssembler(const Element& element,
                                          const FEElementGeometry& feGeometry_,
                                          std::vector<ElementVariables>& elemVars,
                                          const StageParams& stageParams)
    : element_(element)
    , feGeometry_(feGeometry_)
    , elemVariables_(elemVars)
    , stageParams_(stageParams)
    , elementIsGhost_(element.partitionType() == Dune::GhostEntity)
    {}

    /*!
     * \brief Evaluate the complete local residual for the current element.
     */
    ElementResidualVector evalLocalResidual() const
    {
        ElementResidualVector residual(feGeometry_.feBasisLocalView().size());
        residual = 0.0;

        // residual of ghost elements is zero
        if (elementIsGhost_)
            return residual;

        for (std::size_t k = 0; k < stageParams_.size(); ++k)
        {
            LocalOperator localOperator(element_, feGeometry_, elemVariables_[k]);

            if (!stageParams_.skipTemporal(k))
                residual.axpy(stageParams_.temporalWeight(k), localOperator.evalStorage());

            if (!stageParams_.skipSpatial(k))
                residual.axpy(stageParams_.spatialWeight(k), localOperator.evalFluxesAndSources());
        }

        return residual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds
     *        them to the global matrix. The element residual is written into the
     *        right hand side.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(JacobianMatrix& jac,
                                     ResidualVector& res,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        const auto eIdxGlobal = feGeometry_.gridGeometry().elementMapper().index(element_);

        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->evalLocalResidual();
            const auto& localView = feGeometry_.feBasisLocalView();
            for (unsigned int i = 0; i < localView.size(); ++i)
                res[localView.index(i)] += residual[i];
        }
        else if (!elementIsGhost_)
        {
            const auto residual = this->assembleJacobianAndResidual_(jac, partialReassembler);
            const auto& localView = feGeometry_.feBasisLocalView();
            for (unsigned int i = 0; i < localView.size(); ++i)
                res[localView.index(i)] += residual[i];
        }
        else
        {
            const auto& localView = feGeometry_.feBasisLocalView();
            const auto& finiteElement = localView.tree().finiteElement();
            const auto numLocalDofs = finiteElement.localBasis().size();

            for (unsigned int i = 0; i < numLocalDofs; ++i)
            {
                const auto& localKey = finiteElement.localCoefficients().localKey(i);
                if (!isGhostEntity_(localKey.subEntity(), localKey.codim()))
                    continue;

                // set main diagonal entries for the entity
                const auto rowIdx = localView.index(i);

                // TODO: use auto and range based for loop!
                typedef typename JacobianMatrix::block_type BlockType;
                BlockType &J = jac[rowIdx][rowIdx];
                for (int j = 0; j < BlockType::rows; ++j)
                    J[j][j] = 1.0;

                // set residual for the entity
                res[rowIdx] = 0;
            }
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac)
    {
        assembleJacobianAndResidual_(jac);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(ResidualVector& res)
    {
        const auto residual = evalLocalResidual();

        const auto& localView = feGeometry_.feBasisLocalView();
        for (unsigned int i = 0; i < localView.size(); ++i)
            res[localView.index(i)] += residual[i];
    }

protected:

    /*!
     * \brief Returns true if a sub entity of an element is a ghost entity
     */
    bool isGhostEntity_(unsigned int subEntityIdx, unsigned int codim) const
    {
        static constexpr int dim = Element::Geometry::mydimension;

        if (codim == 0)
            return isGhostEntity_(element_.template subEntity<0>(subEntityIdx));
        if (codim == 1)
            return isGhostEntity_(element_.template subEntity<1>(subEntityIdx));
        if constexpr (dim > 1)
            if (codim == 2)
                return isGhostEntity_(element_.template subEntity<2>(subEntityIdx));
        if constexpr (dim > 2)
        {
            assert(codim == 3);
            return isGhostEntity_(element_.template subEntity<3>(subEntityIdx));
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid codimension provided");
    }

    /*!
     * \brief Returns true if an entity is ghost entity
     */
    template<class Entity>
    bool isGhostEntity_(const Entity& e) const
    { return e.partitionType() == Dune::GhostEntity; }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     * \return The element residual at the current solution.
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidual_(JacobianMatrix& A,
                                                       const PartialReassembler* partialReassembler = nullptr)
    {
        using BlockType = typename JacobianMatrix::block_type;
        using PrimaryVariables = typename ElementResidualVector::value_type;
        using Scalar = typename PrimaryVariables::value_type;

        static constexpr unsigned int numPriVars = PrimaryVariables::size();
        static constexpr unsigned int numEq = BlockType::rows;
        static_assert(numPriVars == int(BlockType::cols), "row size mismatch with provided matrix block type.");

        // compute original residuals
        const auto origResiduals = evalLocalResidual();

        // create copy of element solution to undo deflections later
        auto& elemSol = elemVariables_.back().elemSol();
        const auto origElemSol = elemSol;

        ///////////////////////////////////////tsLocalVi///////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        const auto& localView = feGeometry_.feBasisLocalView();
        ElementResidualVector partialDerivs(localView.size());
        for (unsigned int localI = 0; localI < localView.size(); ++localI)
        {
            // dof index and corresponding actual pri vars
            const auto globalI = localView.index(localI);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numPriVars; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the element solution and compute element residual
                    elemSol[localI][pvIdx] = priVar;
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{problem_().paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem_().paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[localI][pvIdx], partialDerivs,
                                                          origResiduals, eps_(elemSol[localI][pvIdx], pvIdx),
                                                          numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (unsigned int localJ = 0; localJ < localView.size(); ++localJ)
                {
                    const auto globalJ = localView.index(localJ);

                    // don't add derivatives for green entities
                    if (!partialReassembler || partialReassembler->dofColor(globalJ) != EntityColor::green)
                    {
                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                        {
                            // A[i][col][eqIdx][pvIdx] is the rate of change of the
                            // the residual of equation 'eqIdx' at dof 'i'
                            // depending on the primary variable 'pvIdx' at dof 'col'
                            A[globalJ][globalI][eqIdx][pvIdx] += partialDerivs[localJ][eqIdx];
                        }
                    }
                }

                // restore the original element solution
                elemSol[localI][pvIdx] = origElemSol[localI][pvIdx];

                // TODO additional dof dependencies
            }
        }

        return origResiduals;
    }

protected:
    //! Return a reference to the underlying problem
    const Problem& problem_() const
    { return elemVariables_.back().gridVariables().problem(); }

private:
    const Element& element_;
    const FEElementGeometry& feGeometry_;
    std::vector<ElementVariables>& elemVariables_;
    const StageParams& stageParams_;
    bool elementIsGhost_;
};

} // end namespace Dumux

#endif
