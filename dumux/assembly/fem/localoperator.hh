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
 * \brief The base class for an element-wise local operator
 *        in finite element schemes.
 */
#ifndef DUMUX_FE_LOCAL_OPERATOR_HH
#define DUMUX_FE_LOCAL_OPERATOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/functions/functionspacebases/subentitydofs.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <dumux/discretization/fem/ipdata.hh>
#include <dumux/discretization/fem/elementboundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise local operator for finite element schemes.
 *        This allows for element-wise evaluation of individual terms
 *        of the equations to be solved.
 * \todo TODO: Doc this!
 * \todo TODO: This currently assumes Standard Galerkin type fem schemes as
 *             well as non-composite function space bases. Can we generalize this?
 */
template<class ElementVariables, class Operators>
class FELocalOperator
{
    // The variables required for the evaluation of the equation
    using GridVars = typename ElementVariables::GridVariables;
    using IpVariables = typename GridVars::IntegrationPointVariables;
    using PrimaryVariables = typename GridVars::PrimaryVariables;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVars::GridGeometry;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    // user-input defined in the problem
    using Problem = typename GridVars::Problem;
    using NumEqVector = typename ProblemTraits<Problem>::NumEqVector;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = NumEqVector::size();

public:
    //! export the grid variables type this residual requires a local view of
    using GridVariables = GridVars;

    //! export a grid-independent alias for compatibility with non grid-based schemes
    using Variables = GridVars;

    //! export underlying scalar type
    using Scalar = typename PrimaryVariables::value_type;

    //! the container storing the residual on all dofs of an element
    using ElementResidualVector = Dune::BlockVector<NumEqVector>;

    //! export a grid-independent alias for compatibility with non grid-based schemes
    using Residual = ElementResidualVector;

    /*!
     * \brief The constructor
     * \note The grid geometry/grid variables local views are expected to
     *       be bound to the same element
     */
    FELocalOperator(const Element& element,
                    const FEElementGeometry& feGeometry,
                    const ElementVariables& elemVars)
    : element_(element)
    , feGeometry_(feGeometry)
    , elemVariables_(elemVars)
    , operators_(element, feGeometry, elemVars)
    {
        setDefaultIntegrationOrders_();
    }

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the terms of the local residual that do not appear in
     *        time derivatives. These are the sources and the fluxes.
     */
    ElementResidualVector evalFluxesAndSources() const
    {
        // evaluate both the volume and the flux terms
        auto volumeTerms = [&] (const auto& ipData, const auto& ipVars)
        { auto s = operators_.source(ipData, ipVars); s *= -1.0; return s; };

        auto fluxTerms = [&] (const auto& ipData, const auto& ipVars)
        { return operators_.flux(ipData, ipVars); };

        auto result = integrateTerms_(volumeTerms, fluxTerms);
        result += evalNeumannSegments_();
        return result;
    }

    /*!
     * \brief Compute the storage term, i.e. the term appearing in the time derivative.
     */
    ElementResidualVector evalStorage() const
    {
        // evaluate both the volume and the flux terms
        auto volumeTerms = [&] (const auto& ipData, const auto& ipVars)
        { return operators_.storage(ipData, ipVars); };

        return integrateVolumeTerms_(volumeTerms);
    }

    ElementResidualVector getEmptyResidual() const
    {
        const auto& localView = feGeometry_.feBasisLocalView();
        ElementResidualVector res(localView.tree().finiteElement().localBasis().size());
        res = 0.0;
        return res;
    }

    // \}

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! \todo TODO: Add interfaces. Or, should this be here at all!?

    //\}

    /*!
     * \name Setter functions for integration orders
     */
    // \{

    /*!
     * \brief Set the integration order to be used for volume integrals.
     * \param intOrder The integration order
     */
    void setIntegrationOrder(unsigned int intOrder)
    {
        intOrder_ = intOrder;
    }

    /*!
     * \brief Set the integration order to be used for boundary integrals.
     * \param intOrderBoundary The boundary integration order
     */
    void setBoundaryIntegrationOrder(unsigned int intOrderBoundary)
    {
        intOrderBoundary_ = intOrderBoundary;
    }

    // \}

protected:
    /*!
     * \brief Integrates the provided terms over the element
     * \todo TODO: Doc requirements on term function interfaces
     */
    template<class VolumeTerms, class FluxTerms>
    ElementResidualVector integrateTerms_(VolumeTerms&& volumeTerms,
                                          FluxTerms&& fluxTerms) const
    {
        const auto& basisLocalView = feGeometry_.feBasisLocalView();
        const auto& localBasis = basisLocalView.tree().finiteElement().localBasis();
        const auto numLocalDofs = localBasis.size();

        ElementResidualVector result(numLocalDofs);
        result = 0.0;

        const auto& geometry = element_.geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder_);
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            FEIntegrationPointData ipData(geometry, quadPoint.position(), localBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            IpVariables ipVars;
            ipVars.update(elemVariables_.elemSol(), problem_(), element_, ipData);

            // evaluate terms and add entries to result
            const auto volume = volumeTerms(ipData, ipVars);
            const auto flux = fluxTerms(ipData, ipVars);

            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < numLocalDofs; ++i)
                    result[i][eqIdx] += qWeight*ipVars.extrusionFactor()
                                        *(volume[eqIdx]*ipData.shapeValue(i)
                                          + flux[eqIdx]*ipData.gradN(i));
        }

        return result;
    }

    /*!
     * \brief Integrates the provided terms over the element
     * \todo TODO: Doc requirements on term function interfaces
     */
    template<class VolumeTerms>
    ElementResidualVector integrateVolumeTerms_(VolumeTerms&& volumeTerms) const
    {
        const auto& basisLocalView = feGeometry_.feBasisLocalView();
        const auto& localBasis = basisLocalView.tree().finiteElement().localBasis();
        const auto numLocalDofs = localBasis.size();

        ElementResidualVector result(numLocalDofs);
        result = 0.0;

        const auto& geometry = element_.geometry();
        const auto& quadRule = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), intOrder_);
        for (const auto& quadPoint : quadRule)
        {
            // Obtain and store shape function values and gradients at the current quad point
            FEIntegrationPointData ipData(geometry, quadPoint.position(), localBasis);

            // calculate secondary variables for the previous and the current solution at the ip
            IpVariables ipVars;
            ipVars.update(elemVariables_.elemSol(), problem_(), element_, ipData);

            // evaluate terms and add entries to result
            const auto volume = volumeTerms(ipData, ipVars);
            Scalar qWeight = quadPoint.weight()*geometry.integrationElement(quadPoint.position());
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                for (unsigned int i = 0; i < numLocalDofs; ++i)
                    result[i][eqIdx] += qWeight*ipVars.extrusionFactor()
                                        *volume[eqIdx]*ipData.shapeValue(i);
        }

        return result;
    }

    /*!
     * \brief Compute the contributions of Neumann boundary conditions.
     */
    ElementResidualVector evalNeumannSegments_() const
    {
        const auto& basisLocalView = feGeometry_.feBasisLocalView();
        const auto& localBasis = basisLocalView.tree().finiteElement().localBasis();
        const auto numLocalDofs = localBasis.size();

        ElementResidualVector result(numLocalDofs);
        result = 0.0;

        // integrate Neumann boundary contribution
        const auto geometry = element_.geometry();
        for (const auto& is : intersections(feGeometry_.gridGeometry().gridView(), element_))
        {
            // only handle faces on the boundary
            if (!is.boundary())
                continue;

            // If none of the dofs living on this intersection has neumann BC defined, skip rest
            // TODO: TIME???
            const auto bcTypes = problem_().boundaryTypes(element_, is);
            if (!bcTypes.hasNeumann())
                continue;

            // get the dofs of this intersection
            const auto isDofs = Dune::Functions::subEntityDOFs(basisLocalView, is);

            // select quadrature rule for intersection faces (dim-1)
            auto insideGeom = is.geometryInInside();
            const auto& faceRule = Dune::QuadratureRules<Scalar, dim-1>::rule(insideGeom.type(), intOrderBoundary_);
            for (const auto& quadPoint : faceRule)
            {
                // position of quadrature point in local coordinates of inside element
                auto local = insideGeom.global(quadPoint.position());

                // evaluate basis functions of all element vertices for quadrature point
                FEIntegrationPointData ipData(geometry, local, localBasis);

                // evaluate primary/secondary variables at integration point
                IpVariables ipVars;
                ipVars.update(elemVariables_.elemSol(), problem_(), element_, ipData);

                // evaluate neumann boundary condition
                const auto neumannFlux = problem_().neumann(element_, is, elemVariables_, ipData, ipVars);

                // get quadrature rule weight for intersection
                Scalar qWeight = quadPoint.weight();
                qWeight *= is.geometry().integrationElement(quadPoint.position());
                qWeight *= ipVars.extrusionFactor();

                // add entries to residual vector
                for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (!bcTypes.isNeumann(eqIdx))
                        continue;

                    for (unsigned int i = 0; i < numLocalDofs; ++i)
                        if (isDofs.contains(i))
                            result[i][eqIdx] += ipData.shapeValue(i)*qWeight*neumannFlux[eqIdx];
                }
            }
        }

        return result;
    }

    //! return reference to the underlying problem
    const Problem& problem_() const
    { return elemVariables_.gridVariables().problem(); }

private:
    //! obtains the integration orders from the input file or default values
    void setDefaultIntegrationOrders_()
    {
        static const auto basisOrder = feGeometry_.feBasisLocalView().tree().finiteElement().localBasis().order();
        static const auto intOrder = getParamFromGroup<unsigned int>(problem_().paramGroup(), "Assembly.FEIntegrationOrder", basisOrder+1);
        static const auto bOrder = getParamFromGroup<unsigned int>(problem_().paramGroup(), "Assembly.FEBoundaryIntegrationOrder", intOrder);

        intOrder_ = intOrder;
        intOrderBoundary_ = bOrder;
    }

    const Element& element_;                 //!< pointer to the element for which the residual is computed
    const FEElementGeometry& feGeometry_;    //!< the local view on the finite element grid geometry
    const ElementVariables& elemVariables_;  //!< the local view on the grid variables

    Operators operators_; //!< evaluates storage/flux operators of the actual equation at integration points
    unsigned int intOrder_;               //!< Integration order used for volume integrals
    unsigned int intOrderBoundary_;       //!< Integration order used for integration of boundary conditions
};

} // end namespace Dumux

#endif
