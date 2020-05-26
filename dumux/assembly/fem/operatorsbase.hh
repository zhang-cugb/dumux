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
 * \brief The base class for operators evaluated at integration points.
 */
#ifndef DUMUX_FE_OPERATORS_BASE_HH
#define DUMUX_FE_OPERATORS_BASE_HH

#include <dune/common/fmatrix.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The base class for operators evaluated at integration points.
 * \todo TODO: Doc this!
 * \todo TODO: This currently assumes Standard Galerkin type fem schemes as
 *             well as non-composite function space bases. Can we generalize this?
 */
template<class GridVarsLocalView>
class FEOperatorsBase
{
    // The variables required for the evaluation of the equation
    using GridVariables = typename GridVarsLocalView::GridVariables;
    using IpVariables = typename GridVariables::IntegrationPointVariables;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename PrimaryVariables::value_type;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVariables::GridGeometry;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    // user-input defined in the problem
    using Problem = typename GridVariables::Problem;
    using NumEqVector = typename ProblemTraits<Problem>::NumEqVector;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = NumEqVector::size();

public:

    // export the types of the flux/source/storage terms
    using FluxTerm = Dune::FieldMatrix<Scalar, numEq, dim>;
    using SourceTerm = NumEqVector;
    using StorageTerm = NumEqVector;

    /*!
     * \brief The constructor
     * \note The grid geometry/grid variables local views are expected to
     *       be bound to the given element
     */
    FEOperatorsBase(const Element& element,
                    const FEElementGeometry& feGeometry,
                    const GridVarsLocalView& gridVarsLocalView)
    : element_(&element)
    , feGeometry_(&feGeometry)
    , gridVarsLocalView_(&gridVarsLocalView)
    {}

    /*!
     * \name Model specific interfaces
     * \note The following method are the model specific implementations of the
     *       operators appearing in the equation and should be overloaded by th
     *       implementations.
     */
    // \{

    /*!
     * \brief Calculate the storage term of the equation
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
     template<class IpData>
     NumEqVector storage(const IpData& ipData,
                         const IpVariables& ipVars) const
    { DUNE_THROW(Dune::NotImplemented, "This model does not implement a storage method!"); }

    /*!
     * \brief Calculate the source term of the equation
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
     template<class IpData>
     NumEqVector source(const IpData& ipData,
                        const IpVariables& ipVars) const
    {
        NumEqVector source(0.0);

        // some references for convenience
        const auto& e = element();
        const auto& fg = feGeometry();
        const auto& gv = gridVariablesLocalView();

        // add contributions from volume flux sources
        source += problem_().source(e, fg, gv, ipData, ipVars);

        // TODO: add contribution from possible point sources

        return source;
    }

    /*!
     * \brief Calculate the flux term of the equation
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
    template<class IpData>
    FluxTerm flux(const IpData& ipData,
                  const IpVariables& ipVars) const
    { DUNE_THROW(Dune::NotImplemented, "This model does not implement a flux method!"); }

    /*!
     * \name Return functions for underlying data structures
     */
    // \{

    /*!
     * \brief Return a reference to the underlying grid element
     */
    const Element& element() const
    { return *element_; }

    /*!
     * \brief Return the local view on the grid geometry
     */
    const FEElementGeometry& feGeometry() const
    { return *feGeometry_; }

    /*!
     * \brief Return reference to the grid variables
     */
    const GridVarsLocalView& gridVariablesLocalView() const
    { return *gridVarsLocalView_; }

    // \}

protected:

    //! return reference to the underlying problem
    const Problem& problem_() const
    { return gridVariablesLocalView().gridVariables().problem(); }

private:

    const Element* element_;                     //!< pointer to the element for which the residual is computed
    const FEElementGeometry* feGeometry_;        //!< the local view on the finite element grid geometry
    const GridVarsLocalView* gridVarsLocalView_; //!< the local view on the grid variables
};

} // end namespace Dumux

#endif
