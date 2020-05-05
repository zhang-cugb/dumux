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
 * \ingroup OneTests
 * \brief The spatial parameters class blood flow problem
 */
#ifndef DUMUX_TEST_1D3D_MICROVESSEL_BlOOD_FLOW_SPATIALPARAMS_HH
#define DUMUX_TEST_1D3D_MICROVESSEL_BlOOD_FLOW_SPATIALPARAMS_HH

#include <vector>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup OneTests
 * \brief Definition of the spatial parameters for the blood flow problem
 */
template<class GridGeometry, class Scalar>
class BloodFlowSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar, BloodFlowSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BloodFlowSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const Scalar r = radius(scv.dofIndex());
        const Scalar gamma = 2.0; // quadratic velocity profile (Poiseuille flow)
        return r*r/(2*(2+gamma));
    }

    /*!
     * \brief Return the radius of the circular pipe for the current sub-control volume in [m].
     */
    Scalar radius(std::size_t eIdxGlobal) const
    { return radius_[eIdxGlobal];}

    /*!
     * \brief Function for defining the velocity estimate.
     */
    Scalar velocityEstimate(std::size_t eIdxGlobal) const
    { return velocityEstimate_[eIdxGlobal]; }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! Set the radii from outside
    void setRadii(const std::vector<Scalar>& radius)
    { radius_ = radius; }

    //! Get the radii for e.g. output
    const std::vector<Scalar>& getRadii() const
    { return radius_; }

    //! Read params from dgf
    template<class GridData>
    void readGridParams(const GridData& gridData)
    {
        const auto& gg = this->gridGeometry();
        auto numElements = gg.gridView().size(0);
        radius_.resize(numElements);
        velocityEstimate_.resize(numElements);

        // gridview is a leafGridView. Parameters are only set on level 0.
        // elements have to inherit spatial parameters from their father.
        for (const auto& element : elements(gg.gridView()))
        {
            auto level0element = element;
            while (level0element.hasFather())
                level0element = level0element.father();

            auto eIdx = gg.elementMapper().index(element);
            radius_[eIdx] = gridData.parameters(level0element)[0];
            velocityEstimate_[eIdx] = gridData.parameters(level0element)[1];
        }
    }

private:
    std::vector<Scalar> radius_;
    std::vector<Scalar> velocityEstimate_;
};

} // end namespace Dumux

#endif
