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
 * \brief The spatial parameters class for the test problem using the 1p cc model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class OnePTestSpatialParams; // TODO rename

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(OnePTestSpatialParams);

}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p cc model
 */
template<class TypeTag>
class OnePTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexSet = typename GridView::IndexSet;
    using ScalarVector = std::vector<Scalar>;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePTestSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {
        permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);

//        lensLowerLeft_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, SpatialParams.LensLowerLeft);
//        lensUpperRight_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, SpatialParams.LensUpperRight);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
//        if (isInLens_(scv.dofPosition()))
//        {
//                return permeabilityLens_;
//        }
//        else
            return permeability_;
    }

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }


private:
//    bool isInLens_(const GlobalPosition &globalPos) const
//    {
//        for (int i = 0; i < dimWorld; ++i) {
//            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
//                return false;
//        }
//        return true;
//    }

//    GlobalPosition lensLowerLeft_;
//    GlobalPosition lensUpperRight_;

    Scalar permeability_;
//    Scalar permeabilityLens_;

    static constexpr Scalar eps_ = 1.0e-7;
};

} // end namespace

#endif
