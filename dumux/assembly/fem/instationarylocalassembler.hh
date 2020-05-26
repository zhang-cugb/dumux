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
template<class Assembler>
class FEInstationaryLocalAssembler
{
    using SolutionVector = typename Assembler::SolutionVector;
    using GridVariables = typename Assembler::GridVariables;
    using GridGeometry = typename GridVariables::GridGeometry;
    using Problem = typename GridVariables::Problem;

    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridVarsLocalView = typename GridVariables::LocalView;

    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using ElementSolution = FEElementSolution<FEElementGeometry, PrimaryVariables>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numEq = PrimaryVariables::size();
    using ModelLocalOperator = typename Assembler::LocalOperator;

public:
    //! TODO: This is contradictory!? Assembler exports local operator, but local assembler exports another one!?
    using LocalOperator = MultiStageLocalOperator<ModelLocalOperator>;
    using ElementResidualVector = typename ModelLocalOperator::ElementResidualVector;

    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ResidualVector = typename Assembler::ResidualVector;

    /*!
     * \brief Construct the local assembler from an assembler.
     */
    explicit FEInstationaryLocalAssembler(const Assembler& assembler)
    : assembler_(assembler)
    {}

    /*!
     * \brief Bind this local assembler to an element of the grid.
     * \param element The grid element
     * \param gridVars An instance of the grid variables
     */
    void bind(const Element& element,
              const GridVariables& gridVars)
    {
        elementIsGhost_ = element.partitionType() == Dune::GhostEntity;

        const auto& prevStageGridVars = assembler().previousStageVariables();
        const auto numStages = prevStageGridVars.size() + 1;

        feGeometry_.clear(); feGeometry_.reserve(numStages);
        gridVarsLocalView_.clear(); gridVarsLocalView_.reserve(numStages);
        modelOperators_.clear(); modelOperators_.reserve(numStages);

        // lambda to add local views for grid vars
        auto addLocalViews = [&] (const auto& curGridVars)
        {
            feGeometry_.emplace_back(curGridVars.gridGeometry());
            gridVarsLocalView_.emplace_back(curGridVars);

            feGeometry_.back().bind(element);
            gridVarsLocalView_.back().bind(element, feGeometry_.back());

            modelOperators_.emplace_back(element, feGeometry_.back(), gridVarsLocalView_.back());
        };

        //! TODO: Do we really need to build geometries for all stages?
        //!       Can we check for adaptivity within a time step in some way?
        //!       Do we need adaptivity during time integration at all or is it
        //!       enough to update after a time step?
        //!       This wouldn't work anyway if the grid changes within stages!?
        //!       Then there is different dofs per stage and it would not be clear
        //!       how to assemble the overall residual for a dof?
        for (const auto& curVars : prevStageGridVars)
            addLocalViews(curVars);
        addLocalViews(gridVars);

        // make multi stage local operator
        // TODO: Pass stage params to bind?! Should assembler have that interface!?
        localOperator_ = std::make_unique<LocalOperator>(modelOperators_, *assembler().stageParams());
    }

    /*!
     * \brief Evaluate the complete local residual for the current element.
     */
    ElementResidualVector evalLocalResidual() const
    {
        // residual of ghost elements is zero
        if (elementIsGhost_)
        {
            ElementResidualVector residual(feGeometry().feBasisLocalView().size());
            residual = 0.0;
            return residual;
        }

        // stationary assembly -> only assemble fluxes and sources
        return localOperator_->evalLocalResidual();
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
        const auto& element = element_();
        const auto eIdxGlobal = feGeometry().gridGeometry().elementMapper().index(element);

        if (partialReassembler && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->evalLocalResidual();
            const auto& localView = feGeometry().feBasisLocalView();
            for (unsigned int i = 0; i < localView.size(); ++i)
                res[localView.index(i)] += residual[i];
        }
        else if (!elementIsGhost_)
        {
            const auto residual = this->assembleJacobianAndResidual_(jac, partialReassembler);
            const auto& localView = feGeometry().feBasisLocalView();
            for (unsigned int i = 0; i < localView.size(); ++i)
                res[localView.index(i)] += residual[i];
        }
        else
        {
            const auto& localView = feGeometry().feBasisLocalView();
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

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            const auto& elemSol = gridVarsLocalView().elemSol();
            res[dofIdx][eqIdx] = elemSol[localDofIdx][pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[dofIdx];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;

            // TODO: Periodic constraints
        };

        enforceDirichletConstraints_(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac)
    {
        assembleJacobianAndResidual_(jac);

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            auto& row = jac[dofIdx];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;
        };

        enforceDirichletConstraints_(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(ResidualVector& res)
    {
        const auto residual = evalLocalResidual();

        const auto& localView = feGeometry().feBasisLocalView();
        for (unsigned int i = 0; i < localView.size(); ++i)
            res[localView.index(i)] += residual[i];

        auto applyDirichlet = [&] (const auto& dirichletValues,
                                   const auto localDofIdx,
                                   const auto dofIdx,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            const auto& elemSol = gridVarsLocalView().elemSol();
            res[dofIdx][eqIdx] = elemSol[localDofIdx][pvIdx] - dirichletValues[pvIdx];
        };

        enforceDirichletConstraints_(applyDirichlet);
    }

    //! The assembler
    const Assembler& assembler() const
    { return assembler_; }

    //! The finite volume geometry
    FEElementGeometry& feGeometry()
    { return feGeometry_.back(); }

    //! The finite volume geometry
    const FEElementGeometry& feGeometry() const
    { return feGeometry_.back(); }

    //! The local view on the grid variables
    GridVarsLocalView& gridVarsLocalView()
    { return gridVarsLocalView_.back(); }

    //! The local view on the grid variables
    const GridVarsLocalView& gridVarsLocalView() const
    { return gridVarsLocalView_.back(); }

    //! Return the underlying local residual
    const LocalOperator& localOperator() const
    { return localOperator_; }

protected:
    //! Return element to which this local assembler is bound
    const Element& element_() const
    { return modelOperators_.back().element(); }

    //! Enforce Dirichlet constraints
    template<typename ApplyFunction>
    void enforceDirichletConstraints_(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        evalDirichletBoundaries_(applyDirichlet);
        // TODO: internal constraints!
        // this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Evaluates Dirichlet boundaries
     */
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries_(ApplyDirichletFunctionType applyDirichlet)
    {
        const auto& elemBcTypes = modelOperators_.back().elemBcTypes();
        // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        if (elemBcTypes.hasDirichlet())
        {
            const auto& localView = feGeometry().feBasisLocalView();
            const auto& finiteElement = localView.tree().finiteElement();

            for (unsigned int localDofIdx = 0; localDofIdx < localView.size(); localDofIdx++)
            {
                const auto& bcTypes = elemBcTypes[localDofIdx];
                if (!bcTypes.hasDirichlet())
                    continue;

                const auto dofIdx = localView.index(localDofIdx);
                const auto& localKey = finiteElement.localCoefficients().localKey(localDofIdx);
                const auto subEntity = localKey.subEntity();
                const auto codim = localKey.codim();

                // values of dirichlet BCs
                PrimaryVariables dirichletValues = getDirichletValues_(subEntity, codim);

                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (bcTypes.isDirichlet(eqIdx))
                    {
                        const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                        assert(0 <= pvIdx && pvIdx < numEq);
                        applyDirichlet(dirichletValues, localDofIdx, dofIdx, eqIdx, pvIdx);
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns the Dirichlet boundary conditions for a sub entity of the element
     */
    PrimaryVariables getDirichletValues_(unsigned int subEntityIdx, unsigned int codim)
    {
        static constexpr int dim = Element::Geometry::mydimension;

        if (codim == 0)
            return getDirichletValues_(element_().template subEntity<0>(subEntityIdx));
        if (codim == 1)
            return getDirichletValues_(element_().template subEntity<1>(subEntityIdx));
        if constexpr (dim > 1)
            if (codim == 2)
                return getDirichletValues_(element_().template subEntity<2>(subEntityIdx));
        if constexpr (dim > 2)
        {
            assert(codim == 3);
            return getDirichletValues_(element_().template subEntity<3>(subEntityIdx));
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid codimension provided");
    }

    /*!
     * \brief Returns the Dirichlet boundary conditions for a sub entity of the element
     */
    template<class SubEntity>
    PrimaryVariables getDirichletValues_(const SubEntity& subEntity)
    {
        return problem_().dirichlet(element_(),
                                    subEntity,
                                    gridVarsLocalView().gridVariables().timeLevel());
    }

    /*!
     * \brief Returns true if a sub entity of an element is a ghost entity
     */
    bool isGhostEntity_(unsigned int subEntityIdx, unsigned int codim) const
    {
        static constexpr int dim = Element::Geometry::mydimension;

        if (codim == 0)
            return isGhostEntity_(element_().template subEntity<0>(subEntityIdx));
        if (codim == 1)
            return isGhostEntity_(element_().template subEntity<1>(subEntityIdx));
        if constexpr (dim > 1)
            if (codim == 2)
                return isGhostEntity_(element_().template subEntity<2>(subEntityIdx));
        if constexpr (dim > 2)
        {
            assert(codim == 3);
            return isGhostEntity_(element_().template subEntity<3>(subEntityIdx));
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
        // TODO: implement this and get rid of if statement
        //       -> first, issue with how diff method enters has to be adressed (see private variable)
        if (diffMethod_ != DiffMethod::numeric)
            DUNE_THROW(Dune::NotImplemented, "analytic differentiation for FEM");

        using BlockType = typename JacobianMatrix::block_type;
        using PrimaryVariables = typename ElementResidualVector::value_type;
        using Scalar = typename PrimaryVariables::value_type;

        static constexpr unsigned int numPriVars = PrimaryVariables::size();
        static constexpr unsigned int numEq = BlockType::rows;
        static_assert(numPriVars == int(BlockType::cols), "row size mismatch with provided matrix block type.");

        // compute original residuals
        const auto origResiduals = evalLocalResidual();

        // create copy of element solution to undo deflections later
        auto& elemSol = gridVarsLocalView().elemSol();
        const auto origElemSol = elemSol;

        ///////////////////////////////////////tsLocalVi///////////////////////////////////////////////////////
        // Calculate derivatives of the residual of all dofs in element with respect to themselves. //
        //////////////////////////////////////////////////////////////////////////////////////////////

        const auto& localView = feGeometry().feBasisLocalView();
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
    { return gridVarsLocalView().gridVariables().problem(); }

private:
    const Assembler& assembler_;          //!< reference to assembler instance

    //! element-local views on the geometry and variables for all stages
    //! these define the state of this class after calling bind()
    std::vector<FEElementGeometry> feGeometry_;
    std::vector<GridVarsLocalView> gridVarsLocalView_;

    bool elementIsGhost_; //!< whether the element's partitionType is ghost
    std::vector<ModelLocalOperator> modelOperators_; //!< TODO: DOC
    std::unique_ptr<LocalOperator> localOperator_; //!< operator evaluating the terms of the equations per element

    DiffMethod diffMethod_{DiffMethod::numeric}; //!< the differentiation method (numeric, analytic, ...)
    // TODO: how should diff method come in? Upon construction or via setter? Or as argument to assemble() ?
};

} // end namespace Dumux

#endif
