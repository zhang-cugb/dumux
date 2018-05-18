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
 * \ingroup MultiDomain
 * \ingroup PoroElastic
 * \ingroup PorousMediumFlow
 * \brief \copydoc Dumux::PMFlowFlowPoroMechanicsCouplingManager
 */

#ifndef DUMUX_PMFLOW_POROMECHANICS_COUPLING_MANAGER_HH
#define DUMUX_PMFLOW_POROMECHANICS_COUPLING_MANAGER_HH

#include <algorithm>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/common/timeloop.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup PoroElastic
 * \ingroup PorousMediumFlow
 * \brief Coupling manager for porous medium flow problems coupled to a poro-mechanical
 *        problem both defined on the same grid. Coupling occurs via the change of porosity
 *        and permeability due to mechanical deformations and the influence of the pore
 *        pressure on the effecive stresses acting on the porous medium.
 *
 * \tparam PMFlowId The domain id of the porous medium flow problem
 * \tparam PoroMechanicalId the domain id of the poro-mechanical problem
 */
template< class MDTraits,
          std::size_t PMFlowId = 0,
          std::size_t PoroMechanicalId = PMFlowId+1 >
class PMFlowFlowPoroMechanicsCouplingManager : public CouplingManager< MDTraits >
{
    using SolutionVector = typename MDTraits::SolutionVector;

    using PMFlowDomainIdType = typename MDTraits::template DomainIdx<PMFlowId>;
    using PoroMechDomainIdType = typename MDTraits::template DomainIdx<PoroMechanicalId>;

    static constexpr auto pmFlowDomainId = PMFlowDomainIdType();
    static constexpr auto poroMechDomainId = PoroMechDomainIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using Scalar = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Scalar);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, NumEqVector);
    template<std::size_t id> using SubDomainSolutionVector = typename GET_PROP_TYPE(SubDomainTypeTag<id>, SolutionVector);
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);
    template<std::size_t id> using ElementBoundaryTypes = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ElementBoundaryTypes);

    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using PrimaryVariables = typename GridVariables<id>::PrimaryVariables;
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache::LocalView;
    template<std::size_t id> using FluxVariablesCache = typename ElementFluxVariablesCache<id>::FluxVariablesCache;
    template<std::size_t id> using FVGridGeometry = typename GridVariables<id>::GridGeometry;
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    // Poro-mechanical problems are instationary by nature
    using TimeLoop = TimeLoopBase< Scalar<pmFlowDomainId> >;

    //! we assume that the two sub-problem operate on the same grid
    static_assert(std::is_same< GridView<pmFlowDomainId>, GridView<poroMechDomainId> >::value,
                  "The grid types of the two sub-problems have to be equal!");

    //! this coupling manager is for cc - box only
    static_assert(FVGridGeometry<poroMechDomainId>::discMethod == DiscretizationMethod::box,
                  "Poro-mechanical problem must be discretized with the box scheme for this coupling manager!");

    static_assert(FVGridGeometry<pmFlowDomainId>::discMethod == DiscretizationMethod::cctpfa ||
                  FVGridGeometry<pmFlowDomainId>::discMethod == DiscretizationMethod::ccmpfa,
                  "Porous medium flow problem must be discretized with a cell-centered scheme for this coupling manager!");

    //! types used for coupling element stencils
    template<std::size_t id>
    using CouplingIndexType = typename std::conditional< id == PMFlowId,
                                                               IndexType<poroMechDomainId>,
                                                               IndexType<pmFlowDomainId> >::type;

    // the coupling context of an element of the poromechanical sub-domain
    struct PoroMechCouplingContext
    {
        // We need unique ptrs because the local views have no default constructor
        Element<pmFlowDomainId> pmFlowElement;
        std::unique_ptr< FVElementGeometry<pmFlowDomainId> > pmFlowFvGeometry;
        std::unique_ptr< ElementVolumeVariables<pmFlowDomainId> > pmFlowElemVolVars;
    };

    // the coupling context of an element of the pm flow sub-domain
    struct PMFlowCouplingContext
    {
        // We need the previous element volume variables of "our own" domain
        // in order to be able to compute the storage terms which form part
        // of the coupling residual
        std::unique_ptr< ElementVolumeVariables<pmFlowDomainId> > pmFlowPrevElemVolVars;
    };

public:

    //! types used for coupling stencils
    //! the poro-mechanical domain elements only couple to themselves on the cc grid
    template<std::size_t id>
    using CouplingStencilType = typename std::conditional< id == PMFlowId,
                                                           std::vector< CouplingIndexType<id> >,
                                                           std::array< CouplingIndexType<id>, 1> >::type;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param pmFlowProblem The porous medium flow problem
     * \param poroMechanicalProblem The poro-mechanical problem
     * \param curSol The current solution
     * \param timeLoop The time loop in case an instationary problem is considered
     */
    void init(std::shared_ptr< Problem<pmFlowDomainId> > pmFlowProblem,
              std::shared_ptr< Problem<poroMechDomainId> > poroMechanicalProblem,
              const SolutionVector& curSol)
    {
        problemTuple_ = std::make_tuple(pmFlowProblem, poroMechanicalProblem);
        curSol_ = curSol;

        initializeCouplingMap_();
    }

    /*!
     * \brief updates the current solution.
     */
    void updateSolution(const SolutionVector& sol)
    {
        curSol_ = sol;
    }

    /*!
     * \brief Return the coupling stencil for a given pm flow domain element
     */
    const CouplingStencilType<pmFlowDomainId>& couplingStencil(const Element<pmFlowDomainId>& element,
                                                               PMFlowDomainIdType,
                                                               PoroMechDomainIdType) const
    { return pmFlowCouplingMap_[ problem<pmFlowDomainId>().fvGridGeometry().elementMapper().index(element) ]; }

    /*!
     * \brief Return the coupling element stencil for a given poro-mechanical domain element
     */
    const CouplingStencilType<poroMechDomainId> couplingStencil(const Element<poroMechDomainId>& element,
                                                                PoroMechDomainIdType,
                                                                PMFlowDomainIdType) const
    {
        const auto eIdx = problem<pmFlowDomainId>().fvGridGeometry().elementMapper().index(element);
        return CouplingStencilType<poroMechDomainId>{ {eIdx} };
    }

    /*!
     * \brief For the assembly of the porous medium flow problem, we need the
     *        element volume variables of the previous solution in order to be able
     *        to evaluate the storage term (coupling enters the porosity).
     */
    template< class Assembler >
    void bindCouplingContext(PMFlowDomainIdType, const Element<pmFlowDomainId>& element, const Assembler& assembler)
    {
        // only do this if the problem is instationary
        if (!assembler.isStationaryProblem())
        {
            // first reset the context
            pmFlowCouplingContext_.pmFlowPrevElemVolVars.reset(nullptr);

            // prepare the previous element volume variables
            auto fvGeometry = localView( problem<pmFlowDomainId>().fvGridGeometry() );
            auto elemVolVars = localView( assembler.gridVariables(pmFlowDomainId).prevGridVolVars() );

            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, assembler.prevSol()[pmFlowDomainId]);

            pmFlowCouplingContext_.pmFlowPrevElemVolVars = std::make_unique< ElementVolumeVariables<pmFlowDomainId> >(elemVolVars);
        }
    }

    /*!
     * \brief For the assembly of the porous medium flow problem, we don't have to prepare anything
     */
    template< class Assembler >
    void bindCouplingContext(PoroMechDomainIdType, const Element<poroMechDomainId>& element, const Assembler& assembler)
    {
        // first reset the context
        poroMechCouplingContext_.pmFlowFvGeometry.reset(nullptr);
        poroMechCouplingContext_.pmFlowElemVolVars.reset(nullptr);

        // prepare the fvGeometry and the element volume variables
        // these quantities will be used later to obtain the effective pressure
        auto fvGeometry = localView( problem<pmFlowDomainId>().fvGridGeometry() );
        auto elemVolVars = localView( assembler.gridVariables(pmFlowDomainId).curGridVolVars() );

        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, curSol_[pmFlowDomainId]);

        poroMechCouplingContext_.pmFlowElement = element;
        poroMechCouplingContext_.pmFlowFvGeometry = std::make_unique< FVElementGeometry<pmFlowDomainId> >(fvGeometry);
        poroMechCouplingContext_.pmFlowElemVolVars = std::make_unique< ElementVolumeVariables<pmFlowDomainId> >(elemVolVars);
    }

    /*!
     * \brief For pmflow -> poro-mechanics, update the entries in the solution vector
     */
    template< class Assembler >
    void updateCouplingContext(PMFlowDomainIdType, PoroMechDomainIdType,
                               IndexType<poroMechDomainId> globalJ, const PrimaryVariables<poroMechDomainId>& priVarsJ, unsigned int pvIdxJ,
                               const Assembler& assembler)
    {
        curSol_[poroMechDomainId][globalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \brief For poro-mechanics -> pmflow, update the element volume variables
     */
    template< class Assembler >
    void updateCouplingContext(PoroMechDomainIdType, PMFlowDomainIdType,
                               IndexType<pmFlowDomainId> globalJ, const PrimaryVariables<pmFlowDomainId>& priVarsJ, unsigned int pvIdxJ,
                               const Assembler& assembler)
    {
        // communicate the deflected pm flow domain primary variable
        curSol_[pmFlowDomainId][globalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        // now, update the coupling context (i.e. elemVolVars)
        const auto& element = poroMechCouplingContext_.pmFlowElement;
        const auto& fvGeometry = *poroMechCouplingContext_.pmFlowFvGeometry;
        poroMechCouplingContext_.pmFlowElemVolVars->bindElement(element, fvGeometry, curSol_[pmFlowDomainId]);
    }

    /*!
     * \brief Update the pm flow volume variables (i.e. porosity) in the coupling
     *        context after the poro-mechanical solution has been deflected.
     */
    template< class Assembler >
    void updateCouplingContext(PoroMechDomainIdType, PoroMechDomainIdType,
                               IndexType<poroMechDomainId> globalJ, const PrimaryVariables<poroMechDomainId>& priVarsJ, unsigned int pvIdxJ,
                               const Assembler& assembler)
    {
        // communicate the deflected displacement
        curSol_[poroMechDomainId][globalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        // now, update the coupling context (i.e. elemVolVars)
        const auto& element = poroMechCouplingContext_.pmFlowElement;
        const auto& fvGeometry = *poroMechCouplingContext_.pmFlowFvGeometry;
        poroMechCouplingContext_.pmFlowElemVolVars->bindElement(element, fvGeometry, curSol_[pmFlowDomainId]);
    }

    /*!
     * \brief We don't have to update anything after the pmflow solution was
     *        deflected during the assembly of a pm flow domain element.
     */
    template< class Assembler >
    void updateCouplingContext(PMFlowDomainIdType, PMFlowDomainIdType,
                               IndexType<pmFlowDomainId> globalJ, const PrimaryVariables<pmFlowDomainId>& priVarsJ, unsigned int pvIdxJ,
                               const Assembler& assembler)
    { /* do nothing here */ }

    /*!
     * \brief Update the pm flow domain volume variables after the coupling context has been updated.
     *        This is necessary because the deformation enters the porosity of the pm flow domain.
     */
    template< class FluxVariablesCacheContainer >
    void updateSelf(PMFlowDomainIdType,
                    const Element<pmFlowDomainId>& element,
                    const FVElementGeometry<pmFlowDomainId>& fvGeometry,
                    GridVolumeVariables<pmFlowDomainId>& gridVolVars,
                    FluxVariablesCacheContainer& fluxVarsCacheContainer)
    {
        // this would require updating all volume variables within the stencil
        DUNE_THROW(Dune::NotImplemented, "Poromechanics framework with grid volume variables caching");
    }

    /*!
     * \brief Update the pm flow domain volume variables after the coupling context has been updated.
     *        This is necessary because the deformation enters the porosity of the pm flow domain.
     */
    template< class FluxVariablesCacheContainer >
    void updateSelf(PMFlowDomainIdType,
                    const Element<pmFlowDomainId>& element,
                    const FVElementGeometry<pmFlowDomainId>& fvGeometry,
                    ElementVolumeVariables<pmFlowDomainId>& elemVolVars,
                    FluxVariablesCacheContainer& fluxVarsCacheContainer)
    {
        elemVolVars.bind(element, fvGeometry, curSol_[pmFlowDomainId]);
    }

    /*!
     * \brief Update the poro-mechanics volume variables after the coupling context has been updated.
     *        This is necessary because the fluid density is stored in the
     */
    template< class FluxVariablesCacheContainer >
    void updateSelf(PoroMechDomainIdType,
                    const Element<poroMechDomainId>& element,
                    const FVElementGeometry<poroMechDomainId>& fvGeometry,
                    GridVolumeVariables<poroMechDomainId>& gridVolVars,
                    FluxVariablesCacheContainer& fluxVarsCacheContainer)
    {
        // this would require updating all volume variables within the stencil
        DUNE_THROW(Dune::NotImplemented, "Poromechanics framework with grid volume variables caching");
    }

    /*!
     * \brief Update the poro-mechanics volume variables after the coupling context has been updated.
     *        This is necessary because the fluid density is stored in the
     */
    template< class FluxVariablesCacheContainer >
    void updateSelf(PoroMechDomainIdType,
                    const Element<poroMechDomainId>& element,
                    const FVElementGeometry<poroMechDomainId>& fvGeometry,
                    ElementVolumeVariables<poroMechDomainId>& elemVolVars,
                    FluxVariablesCacheContainer& fluxVarsCacheContainer)
    {
        elemVolVars.bind(element, fvGeometry, curSol_[poroMechDomainId]);
    }

    /*!
     * \brief Compute the divergence of u for a given element and position.
     *        This is evaluated during porosity computation in the pm flow sub-problem.
     */
    Scalar<pmFlowDomainId> computeDivU(const Element<pmFlowDomainId>& element,
                                       const typename Element<pmFlowDomainId>::Geometry::GlobalCoordinate& globalPos) const
    {
        const auto& poroMechFvGridGeometry = problem<poroMechDomainId>().fvGridGeometry();
        const auto poroMechElemSol = elementSolution( element, curSol_[poroMechDomainId],  poroMechFvGridGeometry );
        const auto gradU = evalGradients(element, element.geometry(), poroMechFvGridGeometry, poroMechElemSol, globalPos);

        Scalar<pmFlowDomainId> divU = 0.0;
        for (int i = 0; i < GridView<poroMechDomainId>::dimension; ++i)
            divU += gradU[i][i];
        return divU;
    }

    /*!
     * \brief Evaluates the coupling residual of the porous medium flow domain with
     *        the poro-mechanical domain. The deformation might has an effect on both
     *        the permeability as well as the porosity. Thus, we need to compute fluxes
     *        and the storage term for the coupling residual
     */
    typename LocalResidual<pmFlowDomainId>::ElementResidualVector
    evalCouplingResidual(PMFlowDomainIdType,
                         const Element<pmFlowDomainId>& element,
                         const FVElementGeometry<pmFlowDomainId>& fvGeometry,
                         const ElementVolumeVariables<pmFlowDomainId>& elemVolVars,
                         const ElementBoundaryTypes<pmFlowDomainId>& elemBcTypes,
                         const ElementFluxVariablesCache<pmFlowDomainId>& elemFluxVarsCache,
                         LocalResidual<pmFlowDomainId> localResidual,
                         PoroMechDomainIdType,
                         IndexType<poroMechDomainId> globalJ)
    {
        auto res = localResidual.evalFluxAndSource(element, fvGeometry, elemVolVars, elemFluxVarsCache, elemBcTypes);

        // If the residual instationary, evaluate storage
        if (!localResidual.isStationary())
            res += localResidual.evalStorage(element, fvGeometry, *pmFlowCouplingContext_.pmFlowPrevElemVolVars, elemVolVars);

        return res;
    }

    /*!
     * \brief Evaluates the coupling residual of the poromechanical domain with
     *        the porous medium flow domain. The pressure has an effect on the
     *        mechanical stresses as well as the body forces. Thus, we have to
     *        compute the fluxes and the source term here.
     */
    typename LocalResidual<poroMechDomainId>::ElementResidualVector
    evalCouplingResidual(PoroMechDomainIdType,
                         const Element<poroMechDomainId>& element,
                         const FVElementGeometry<poroMechDomainId>& fvGeometry,
                         const ElementVolumeVariables<poroMechDomainId>& elemVolVars,
                         const ElementBoundaryTypes<poroMechDomainId>& elemBcTypes,
                         const ElementFluxVariablesCache<poroMechDomainId>& elemFluxVarsCache,
                         LocalResidual<poroMechDomainId> localResidual,
                         PMFlowDomainIdType,
                         IndexType<pmFlowDomainId> globalJ)
    {
        return localResidual.evalFluxAndSource(element, fvGeometry, elemVolVars, elemFluxVarsCache, elemBcTypes);
    }

    //! Return a const reference to one of the problems
    template<std::size_t id> const Problem<id>& problem() const
    { return *std::get<id>(problemTuple_); }

    //! Return reference to one of the problems
    template<std::size_t id> Problem<id>& problem()
    { return *std::get<id>(problemTuple_); }

    //! Return the coupling context (used in mechanical sub-problem to compute effective pressure)
    const PoroMechCouplingContext& poroMechanicsCouplingContext() const
    { return poroMechCouplingContext_; }

private:
    /*!
     * \brief Initializes the pm flow domain coupling map. Since the elements
     *        of the poro-mechanical domain only couple to the same elements, we
     *        don't set up the maps here but return copy of the "stencil" always.
     */
    void initializeCouplingMap_()
    {
        // some references for convenience
        const auto& pmFlowGridGeometry = problem<pmFlowDomainId>().fvGridGeometry();
        const auto& poroMechGridGeometry = problem<poroMechDomainId>().fvGridGeometry();
        const auto& vertexMapper = pmFlowGridGeometry.vertexMapper();

        // make sure the two grids are really the same. Note that if the two grids
        // happen to have equal number of elements by chance, we don't detect this source of error.
        if (pmFlowGridGeometry.gridView().size(0) != poroMechGridGeometry.gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "The two sub-problems are assumed to work on the same mesh!");

        // set up the coupling stencil
        static constexpr int dim = GridView<pmFlowDomainId>::dimension;
        pmFlowCouplingMap_.resize(pmFlowGridGeometry.gridView().size(0));
        for (const auto& element : elements(pmFlowGridGeometry.gridView()))
        {
            const auto eIdx = pmFlowGridGeometry.elementMapper().index(element);

            // firstly, the element couples to the nodal dofs in itself
            for (int i = 0; i < element.geometry().corners(); ++i)
                pmFlowCouplingMap_[eIdx].push_back( vertexMapper.subIndex(element, i , dim) );

            // the pm flow problem couples to the same elements as in its own stencil
            // due to the dependency of the residual on all permeabilities in its stencil,
            // which in turn depend on the mechanical deformations
            const auto& inverseConnectivity = pmFlowGridGeometry.connectivityMap()[eIdx];
            for (const auto& dataJ : inverseConnectivity)
            {
                const auto elemJ = pmFlowGridGeometry.element(dataJ.globalJ);
                for (int i = 0; i < elemJ.geometry().corners(); ++i)
                    pmFlowCouplingMap_[dataJ.globalJ].push_back( vertexMapper.subIndex(elemJ, i , dim) );
            }
        }

        // make stencils unique
        for (auto& stencil : pmFlowCouplingMap_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
    }

    // tuple for storing pointers to the sub-problems
    using PMFlowProblemPtr = std::shared_ptr< Problem<pmFlowDomainId> >;
    using PoroMechanicalProblemPtr = std::shared_ptr< Problem<poroMechDomainId> >;
    std::tuple<PMFlowProblemPtr, PoroMechanicalProblemPtr> problemTuple_;

    // Container for storing the coupling element stencils for the pm flow domain
    std::vector<CouplingStencilType<pmFlowDomainId>> pmFlowCouplingMap_;

    //! Store copy of the current solution
    SolutionVector curSol_;

    // the coupling contexts
    PMFlowCouplingContext pmFlowCouplingContext_;
    PoroMechCouplingContext poroMechCouplingContext_;
};

} //end namespace Dumux

#endif
