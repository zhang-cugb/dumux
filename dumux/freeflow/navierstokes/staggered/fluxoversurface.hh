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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::FluxOverSurface
 */
#ifndef DUMUX_FLUX_OVER_SURFACE_STAGGERED_HH
#define DUMUX_FLUX_OVER_SURFACE_STAGGERED_HH

#include <numeric>
#include <functional>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/geometry/makegeometry.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief  Class used to calculate fluxes over surfaces. This only works for the staggered grid discretization.
 */
template <class TypeTag>
class FluxOverSurface
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using Element = typename GridView::template Codim<0>::Entity;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto surfaceDim = dimWorld - 1;
    using SurfaceGeometryType = Dune::MultiLinearGeometry< Scalar, surfaceDim, dimWorld >;

    /*!
     * \brief Auxiliary class that holds the surface-specific data.
     */
    template<int mydim, int coordDim, class Scalar=double>
    class SurfaceData
    {
        using Geo = SurfaceGeometryType;
    public:

        SurfaceData() {}

        SurfaceData(Geo&& geo)
        {
            values_.resize(1);
            geometries_.push_back(std::move(geo));
        }

        SurfaceData(std::vector<Geo>&& geo)
        {
            values_.resize(geo.size());
            geometries_ = std::forward<decltype(geo)>(geo);
        }

        void addValue(int surfaceIdx, const CellCenterPrimaryVariables& value)
        {
            values_[surfaceIdx] += value;
        }

        void addSubSurface(Geo&& geo)
        {
            values_.emplace_back(0.0);
            geometries_.push_back(std::move(geo));
        }

        auto& subSurfaces() const
        {
            return geometries_;
        }

        CellCenterPrimaryVariables& value(int surfaceIdx) const
        {
            return values_[surfaceIdx];
        }

        auto& values() const
        {
            return values_;
        }

        void printSurfaceBoundaries(int surfaceIdx) const
        {
            const auto& geometry = geometries_[surfaceIdx];
            for(int i = 0; i < geometry.corners(); ++i)
                std::cout << "(" << geometry.corner(i) << ")  ";
        }

        void resetValues()
        {
            std::fill(values_.begin(), values_.end(), CellCenterPrimaryVariables(0.0));
        }

    private:

        std::vector<Geo> geometries_;
        std::vector<CellCenterPrimaryVariables> values_;
    };


public:

    using SurfaceList = std::vector<SurfaceGeometryType>;

    /*!
     * \brief The constructor
     *
     * \param assembler The assembler
     * \param sol The solution vector
     */
    FluxOverSurface(const Problem& problem,
                    const GridVariables& gridVariables,
                    const SolutionVector& sol)
    : problem_(problem),
      gridVariables_(gridVariables),
      sol_(sol)
    {
        verbose_  = getParamFromGroup<bool>(problem.paramGroup(), "FluxOverSurface.Verbose", false);
    }

    /*!
     * \brief Add a collection of sub surfaces under a given name
     *
     * \param name The name of the surface
     * \param surfaces The list of sub surfaces
     */
    void addSurface(const std::string& name, SurfaceList&& surfaces)
    {
        surfaces_[name] = SurfaceData<surfaceDim,dim>(std::forward<decltype(surfaces)>(surfaces));
    }

    /*!
     * \brief Add a surface under a given name, specifying the surface's corner points.
     *        This is a specialization for 2D, therefore the surface is actually a line.
     *
     * \param name The name of the surface
     * \param p0 The first corner
     * \param p1 The second corner
     */
    void addSurface(const std::string& name, const GlobalPosition& p0, const GlobalPosition& p1)
    {
        surfaces_[name].addSubSurface(makeSurface(std::vector<GlobalPosition>{p0, p1}));
    }

    /*!
     * \brief Add a surface under a given name, specifying the surface's corner points.
     *        This is a specialization for 3D.
     *
     * \param name The name of the surface
     * \param p0 The first corner
     * \param p1 The second corner
     * \param p2 The third corner
     * \param p3 The fourth corner
     */
    void addSurface(const std::string& name,
                    const GlobalPosition& p0,
                    const GlobalPosition& p1,
                    const GlobalPosition& p2,
                    const GlobalPosition& p3)
    {
        surfaces_[name].addSubSurface(makeSurface(std::vector<GlobalPosition>{p0, p1, p2, p3}));
    }

    /*!
     * \brief Creates a geometrical surface object for (2D).
     *
     * \param corners The vector storing the surface's corners
     */
    static SurfaceGeometryType makeSurface(const std::vector<Dune::FieldVector<Scalar, 2>>& corners)
    {
            return SurfaceGeometryType(Dune::GeometryTypes::line, corners);
    }

    /*!
     * \brief Creates a geometrical surface object for (3D).
     *
     * \param corners The vector storing the surface's corners
     */
    static SurfaceGeometryType makeSurface(const std::vector<Dune::FieldVector<Scalar, 3>>& corners)
    {
        return makeDuneQuadrilaterial(corners);
    }

    /*!
     * \brief Calculate the mass or mole fluxes over all surfaces
     */
    void calculateMassOrMoleFluxes()
    {
        auto fluxType = [this](const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& elemFaceVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            LocalResidual localResidual(&problem_);
            return localResidual.computeFluxForCellCenter(problem_, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);
        };

        calculateFluxes(fluxType);
    }

    /*!
     * \brief Calculate the volume fluxes over all surfaces.
     */
    void calculateVolumeFluxes()
    {
        const auto isCompositional = std::integral_constant<bool, (GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents() > 1) >();
        calculateVolumeFluxesImpl_(isCompositional);
    }

    /*!
     * \brief Calculate the fluxes over all surfaces for a given flux type.
     *
     * \param fluxType The flux type. This can be a lambda of the following form:
     *                 [](const auto& element,
                          const auto& fvGeometry,
                          const auto& elemVolVars,
                          const auto& elemFaceVars,
                          const auto& scvf,
                          const auto& elemFluxVarsCache)
                          { return ... ; }
     */
    template<class FluxType>
    void calculateFluxes(const FluxType& fluxType)
    {
        // make sure to reset all the values of the surfaces, in case this method has been called already before
        for(auto&& surface : surfaces_)
            surface.second.resetValues();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(problem_.fvGridGeometry().numFaceDofs(), false);

        auto elemVolVars = localView(gridVariables_.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
        auto elemFaceVars = localView(gridVariables_.curGridFaceVars());

        for(auto&& element : elements(problem_.fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(problem_.fvGridGeometry());
            fvGeometry.bindElement(element);

            elemVolVars.bind(element, fvGeometry, sol_);
            elemFaceVars.bind(element, fvGeometry, sol_);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for(auto && scvf : scvfs(fvGeometry))
            {
                const auto dofIdx = scvf.dofIndex();
                // do nothing of the dof was already visited
                if(dofVisited[dofIdx])
                    continue;

                dofVisited[dofIdx] = true;

                // iterate through all surfaces and check if the flux at the given position
                // should be accounted for in the respective surface
                for(auto&& surface : surfaces_)
                {
                    const auto& subSurfaces = surface.second.subSurfaces();

                    for(int surfaceIdx = 0; surfaceIdx < subSurfaces.size(); ++surfaceIdx)
                    {
                        if(intersectsPointGeometry(scvf.center(), subSurfaces[surfaceIdx]))
                        {
                            const auto result = fluxType(element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);

                            surface.second.addValue(surfaceIdx, result);

                            if(verbose_)
                            {
                                std::cout << "Flux at face "  << scvf.center() << " (" << surface.first << "): " << result;
                                std::cout << " (directionIndex: " << scvf.directionIndex() << "; surface boundaries: ";
                                surface.second.printSurfaceBoundaries(surfaceIdx);
                                std::cout << ", surfaceIdx " << surfaceIdx << ")" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Return the fluxes of the individual sub surface of a given name.
     *
     * \param name The name of the surface
     */
    auto& values(const std::string& name) const
    {
        return surfaces_.at(name).values();
    }

    /*!
     * \brief Return the cumulative net fluxes of a surface of a given name.
     *
     * \param name The name of the surface
     */
    auto netFlux(const std::string& name) const
    {
        const auto& surfaceResults = values(name);
        return std::accumulate(surfaceResults.begin(), surfaceResults.end(), CellCenterPrimaryVariables(0.0));
    }

private:

    /*!
     * \brief Calculate the volume fluxes over all surfaces for compositional models.
     *        This method simply averages the densities between two adjacent cells.
     */
    void calculateVolumeFluxesImpl_(std::true_type)
    {
        auto fluxType = [this](const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& elemFaceVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            LocalResidual localResidual(&problem_);
            const auto massOrMoleFlux = localResidual.computeFluxForCellCenter(problem_, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);

            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
            const auto density = [useMoles](const auto& volVars)
            {
                return useMoles ? volVars.molarDensity() : volVars.density() ;
            };

            const auto avgDensity = 0.5*density(insideVolVars) + 0.5*density(outsideVolVars);

            constexpr auto replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
            constexpr auto numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents();

            const Scalar cumulativeFlux = [replaceCompEqIdx, numComponents, &massOrMoleFlux]()
            {
                Scalar result = 0.0;
                if(replaceCompEqIdx < numComponents)
                    result = massOrMoleFlux[replaceCompEqIdx];
                else
                {
                    for(int i = 0; i < numComponents; ++i)
                        result += massOrMoleFlux[i];
                }

                return result;
            }();

            CellCenterPrimaryVariables tmp(0.0);
            tmp[0] = cumulativeFlux / avgDensity;
            return tmp;
        };

        calculateFluxes(fluxType);
    }

    /*!
     * \brief Calculate the volume fluxes over all surfaces for non-compositional models.
     *        This method simply averages the densities between two adjacent cells.
     */
    void calculateVolumeFluxesImpl_(std::false_type)
    {
        auto fluxType = [this](const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& elemFaceVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            LocalResidual localResidual(&problem_);
            const Scalar totalMassFlux = localResidual.computeFluxForCellCenter(problem_, element, fvGeometry, elemVolVars,
                                                                                elemFaceVars, scvf, elemFluxVarsCache)[0];

            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            const auto avgDensity = 0.5*insideVolVars.density() + 0.5*outsideVolVars.density();

            CellCenterPrimaryVariables tmp(0.0);
            tmp[0] = totalMassFlux / avgDensity;
            return tmp;
        };

        calculateFluxes(fluxType);
    }

    std::map<std::string, SurfaceData<surfaceDim ,dim> > surfaces_;
    const Problem& problem_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    bool verbose_;
};

} //end namespace

#endif