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
/**
 * @file
 * @brief  Definition of a simple Stokes problem
 */
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>
#include <typeinfo>

namespace Dumux
{

template <class TypeTag>
class StokesTestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
// SET_TYPE_PROP(StokesTestProblem, GridCreator, MatchingGridCreator<TypeTag, 3>);
#else
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
// SET_TYPE_PROP(StokesTestProblem, GridCreator, MatchingGridCreator<TypeTag, 2>);
#endif

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, StokesTestProblem<TypeTag>);

SET_PROP(StokesTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> > type;
};


SET_BOOL_PROP(StokesTestProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(StokesTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(StokesTestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(StokesTestProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(StokesTestProblem, EnableInertiaTerms, false);

// Set the grid parameter group
SET_STRING_PROP(StokesTestProblem, GridParameterGroup, "StokesGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);

}

/*!
 * \ingroup StaggeredStokesModel
 * \ingroup ImplicitTestProblems
 * \brief Stokes flow problem with nitrogen (N2) flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are set to Dirichlet with outflow on the right boundary. The mass balance has
 * outflow bcs, which are replaced in the localresidual by the sum
 * of the momentum balance equations in case of Dirichlet bcs for the momentum balance.
 * In the middle of the right boundary, one vertex receives Dirichlet bcs to set the pressure level.
 * The flow velocity starts with 0 m/s. A flow field evolves with a maximum velocity, which is
 * varied time-dependently using a sinus function and a period of 3000s.
 *
 * This problem uses the \ref StaggeredStokesModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes</tt>
 */
template <class TypeTag>
class StokesTestProblem : public NavierStokesProblem<TypeTag>
{
    typedef NavierStokesProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag) GlobalTypeTag;
    typedef typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager) CouplingManager;
    typedef typename GET_PROP_TYPE(GlobalTypeTag, StokesProblemTypeTag) StokesProblemTypeTag;

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

public:
    StokesTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        vIn_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Velocity);

        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, StokesGrid, LowerLeft)[0];
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, StokesGrid, UperRight)[0];

        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, StokesGrid, LowerLeft)[1];
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, StokesGrid, UpperRight)[1];
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {
        return name_+"_stokes";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return 273.15 + 10.00; // -> 10C
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

//        const Scalar time = this->time.Manager().time();

        // Fetzer2017a
        values.setAllDirichlet();

//        if(onUpperBoundary_(globalPos))
//        	values.setNeumann(transportEqIdx);

        // Left inflow boundaries should be Neumann, otherwise the
        // evaporative fluxes are much more grid dependent
        if(onLeftBoundary_(globalPos))
        {
//        	values.setNeumann(transportEqIdx);

        if(onUpperBoundary_(globalPos)) // corner point
        values.setAllDirichlet();
        }

        if(onRightBoundary_(globalPos))
        {
        values.setAllOutflow();
        if(onUpperBoundary_(globalPos)) // corner point
        values.setAllDirichlet();
        }

        if(onLowerBoundary_(globalPos))
        {
//        	values.setAllDirichlet();
//        	values.setNeumann(transportEqIdx);
//
//        	if(globalPos[0] > runUpDistanceX - eps_ && time > initializationTime_)
//        	{
        values.setAllCouplingDirichlet();
        values.setCouplingNeumann(momentumXBalanceIdx);
        values.setCouplingNeumann(momentumYBalanceIdx);
//        	}
        }

        // the mass balance has to be of type outflow
        // it does not get a coupling condition, since pn is a condition for Stokes
        values.setOutflow(massBalanceIdx);

        // set pressure at one point, do NOT specif this if the Darcy domain
        // has a Dirichlet condition for pressure
        if(onRightBoundary_(globalPos))
        {
//        	if(time > initializationTime_)
//        		values.setDirichlet(pressureIdx);
//        	else
        if(!onLowerBoundary_(globalPos) && !onUpperBoundary_(globalPos))
        values.setDirichlet(pressureIdx);
        }

        // ??
//        // set Dirichlet values for the velocity everywhere
//        values.setDirichlet(momentumBalanceIdx);
//
//        // set a fixed pressure in one cell
//        if (isOutlet(globalPos))
//        {
//            values.setDirichlet(massBalanceIdx);
//            values.setOutflow(momentumBalanceIdx);
//        }
//        else
//            values.setOutflow(massBalanceIdx);

        return values;
    }


//    /*!
//     * \brief Evaluate the boundary conditions for a dirichlet
//     *        control volume.
//     *
//     * \param element The element
//     * \param scvf The subcontrolvolume face
//     */
//    BoundaryValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
//    {
//        const auto& globalPos = scvf.dofPosition();
//        BoundaryValues values;
//        values[pressureIdx] = 0.0;
//        values[velocityXIdx] = 0.0;
//        values[velocityYIdx] = 0.0;
//
//        if(isInlet(globalPos))
//        {
//            values[velocityXIdx] = vIn_;
//            values[velocityYIdx] = 0.0;
//        }
//
//        if(couplingManager().isStokesCouplingEntity(element, scvf))
//        {
//            values[velocityYIdx] = couplingManager().darcyData().boundaryVelocity(scvf); // TODO
//        }
//
//        return values;
//
//    }
    // Fetzer2017a
    BoundaryValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
    BoundaryValues values;

    values[velocityXIdx] = vIn_;
    values[velocityYIdx] = 0.0;
    values[pressureIdx] = 1.0e5;
//    	values[massOrMoleFracIdx] = refMassfrac();

    return values;
    }

    // Fetzer2017a
    BoundaryValues neumannAtPos(const GlobalPosition &globalPos) const
    {
    BoundaryValues values(0.0);

//    	const Scalar xVelocity_(globalPos);
//
//    	if(onLeftBoundary_(globalPos) && gobalPos[1] > bBoxMin_[1] - eps_ && globalPos[1] < bBoxMax_[1] + eps_)
//    	{
//    		// rho*v*X at inflow
//    		values[transportEqIdx] = -1.0 * xVelocity * density * massFracDirichlet_;
//    	}
    return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values;
        values[pressureIdx] = 0.0;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        return values;
    }

    // \}

    /*!
     * \brief Set the coupling manager
     * \param couplingManager The coupling manager
     *
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     *
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:

//    bool isInlet(const GlobalPosition& globalPos) const
//    {
//        return globalPos[0] < this->bBoxMin()[0] + eps_;
//    }
//
//    bool isOutlet(const GlobalPosition& globalPos) const
//    {
//        return globalPos[0] > this->bBoxMax()[0] - eps_;
//    }
//
//    bool isWall(const GlobalPosition& globalPos) const
//    {
//        return globalPos[0] > this->bBoxMin()[0] + eps_ || globalPos[0] < this->bBoxMax()[0] - eps_;
//    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    Scalar eps_;
    Scalar vIn_;
    std::string name_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif
