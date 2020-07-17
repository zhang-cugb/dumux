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
 * \ingroup BoxDiscretization
 * \brief Calculates the element-wise residual for the box scheme
 */
#ifndef DUMUX_FACECENTERED_LOCAL_RESIDUAL_HH
#define DUMUX_FACECENTERED_LOCAL_RESIDUAL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/assembly/fvlocalresidual.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief The element-wise residual for the box scheme
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class FaceCenteredLocalResidual : public FVLocalResidual<TypeTag>
{
    using ParentType = FVLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! evaluate flux residuals for one sub control volume face and add to residual
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const
    {
        const auto flux = evalFlux(problem, element, fvGeometry, elemVolVars, elemBcTypes, elemFluxVarsCache, scvf);
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        residual[insideScv.localDofIndex()] += flux;
    }

    //! evaluate flux residuals for one sub control volume face
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementBoundaryTypes& elemBcTypes,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);

        // inner faces
        if (!scvf.boundary())
        {
            flux += this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        }

        // boundary faces
        else
        {
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

            if (scv.boundary())
            {}
            else if (scvf.isLateral())
            {
                flux += this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            }
        //     const auto& bcTypes = elemBcTypes[scv.localDofIndex()];

        //     // Treat Neumann and Robin ("solution dependent Neumann") boundary conditions.
        //     // For Dirichlet there is no addition to the residual here but they
        //     // are enforced strongly by replacing the residual entry afterwards.
        //     if (bcTypes.hasNeumann())
        //     {
        //         auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

        //         // multiply neumann fluxes with the area and the extrusion factor
        //         neumannFluxes *= scvf.area()*elemVolVars[scv].extrusionFactor();

        //         // only add fluxes to equations for which Neumann is set
        //         for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
        //             if (bcTypes.isNeumann(eqIdx))
        //                 flux[eqIdx] += neumannFluxes[eqIdx];
        //     }
        }

        return flux;
    }
};

} // end namespace Dumux

#endif
