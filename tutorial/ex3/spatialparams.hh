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
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
#ifndef DUMUX_EXERCISE_THREE_SPATIAL_PARAMS_HH
#define DUMUX_EXERCISE_THREE_SPATIAL_PARAMS_HH

// include parent spatialparameters
#include <dumux/material/spatialparams/implicit.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux {
//forward declaration
template<class TypeTag>
class ExerciseThreeSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ExerciseThreeSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ExerciseThreeSpatialParams, SpatialParams,
        ExerciseThreeSpatialParams<TypeTag>);

// Set the material law
SET_PROP(ExerciseThreeSpatialParams, MaterialLaw)
{
private:
    // material law typedefs
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // select material law to be used
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    // adapter for absolute law
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
template<class TypeTag>
class ExerciseThreeSpatialParams: public ImplicitSpatialParams<TypeTag>
{
    // Get informations for current implementation via property system
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum
    {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    // Get object types for function arguments
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

public:
    // get material law from property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    // determine appropriate parameters depending on selected materialLaw
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    const Dune::FieldMatrix<Scalar, dim, dim> &intrinsicPermeability(const Element &element,
                                                    const FVElementGeometry &fvGeometry,
                                                    const int scvIdx) const
    {
        if (isInLens(element.geometry().center()))
            return KLens_;
        return K_;
    }

    /*! Defines the porosity \f$[-]\f$ of the porous medium depending
     * on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function porosityAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        if (isInLens(element.geometry().center()))
            return 0.1;
        return 0.2;
    }

    /*! Returns the parameter object for the material law (i.e. Brooks-Corey)
     *  depending on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function materialLawParamsAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
        if (isInLens(element.geometry().center()))
            return materialParamsLens_;
        return materialParams_;
    }

    // constructor
    ExerciseThreeSpatialParams(const GridView& gridView) :
        ImplicitSpatialParams<TypeTag>(gridView),
        K_(0),
        KLens_(0)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        for (int i = 0; i < dim; i++)
        {
            K_[i][i] = 1e-7;
            KLens_[i][i] = 1e-10;
        }

        //set residual saturations
        materialParams_.setSwr(0.0);
        materialParamsLens_.setSwr(0.1);
        materialParams_.setSnr(0.0);
        materialParamsLens_.setSnr(0.1);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(500.0);
        materialParamsLens_.setPe(1000.0);
        materialParams_.setLambda(2);
        materialParamsLens_.setLambda(2);
    }

    bool isInLens(const Dune::FieldVector<Scalar, dimWorld>& globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return (x < 40 && x > 20 && y > 35 && y < 45) ||
               (x < 50 && x > 30 && y < 30 && y > 15);
    }

private:

    Dune::FieldMatrix<Scalar, dim, dim> K_;
    Dune::FieldMatrix<Scalar, dim, dim> KLens_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;
    MaterialLawParams materialParamsLens_;
};
} // end namespace
#endif