// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief The spatial parameters class for the matrix problem
 */
#ifndef DUMUX_2P2C_MATRIX_SPATIALPARAMS_HH
#define DUMUX_2P2C_MATRIX_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p2c/implicit/properties.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPTwoCMatrixSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPTwoCMatrixSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPTwoCMatrixSpatialParams, SpatialParams, TwoPTwoCMatrixSpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(TwoPTwoCMatrixSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}

template<class TypeTag>
class TwoPTwoCMatrixSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = Scalar;

    TwoPTwoCMatrixSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        K_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MatrixPermeability);

        // residual saturations
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        materialLawParams_.setPe(5e4);
        materialLawParams_.setLambda(2.0);
    }

    /*!
     * \brief Return the intrinsic permeability for the given position in [m^2].
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    { return K_; }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The solution for all dofs of the element
     */
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return 0; }

    /*!
     * \brief Define the porosity in [-].
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    /*!
     * \brief Returns the parameter object for the capillary-pressure/
     *        saturation material law
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialLawParams_; }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    Scalar solidHeatCapacity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolutionVector& elemSol) const
    { return 790; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return 2700; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    { return 2.8; }

private:
    Scalar K_;
    MaterialLawParams materialLawParams_;
};
} //end namespace

#endif
