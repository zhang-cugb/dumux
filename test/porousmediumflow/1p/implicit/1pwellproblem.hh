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
 * \brief A test problem for the one-phase model:
 * A well extracts water with a given pumping rate.
 */
#ifndef DUMUX_1PDUPUITTHIEM_PROBLEM_HH
#define DUMUX_1PDUPUITTHIEM_PROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/linear/amgbackend.hh>

#include "1ptestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePDupuitThiemProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePDupuitThiemProblem, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePTestBoxProblem, INHERITS_FROM(BoxModel, OnePDupuitThiemProblem));
NEW_TYPE_TAG(OnePTestCCProblem, INHERITS_FROM(CCTpfaModel, OnePDupuitThiemProblem));

SET_PROP(OnePDupuitThiemProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
// SET_TYPE_PROP(OnePDupuitThiemProblem, Grid, Dune::YaspGrid<2>);
SET_TYPE_PROP(OnePDupuitThiemProblem, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2> >);
//SET_TYPE_PROP(OnePDupuitThiemProblem, Grid, Dune::UGGrid<2>);
//SET_TYPE_PROP(OnePDupuitThiemProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePDupuitThiemProblem, Problem, Dumux::OnePDupuitThiemProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePDupuitThiemProblem, SpatialParams, Dumux::OnePTestSpatialParams<TypeTag> );

// Linear solver settings
SET_TYPE_PROP(OnePDupuitThiemProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );

NEW_TYPE_TAG(OnePTestBoxProblemWithAMG, INHERITS_FROM(OnePTestBoxProblem));
NEW_TYPE_TAG(OnePTestCCProblemWithAMG, INHERITS_FROM(OnePTestCCProblem));
// Solver settings for the tests using AMG
SET_TYPE_PROP(OnePTestBoxProblemWithAMG, LinearSolver, Dumux::AMGBackend<TypeTag> );
SET_TYPE_PROP(OnePTestCCProblemWithAMG, LinearSolver, Dumux::AMGBackend<TypeTag> );

// if FoamGrid is available, test for dim < dimWorld
#if HAVE_DUNE_FOAMGRID
NEW_TYPE_TAG(OnePOneDThreeDTestProblem, INHERITS_FROM(OnePDupuitThiemProblem));
NEW_TYPE_TAG(OnePOneDThreeDTestBoxProblem, INHERITS_FROM(BoxModel, OnePOneDThreeDTestProblem));
NEW_TYPE_TAG(OnePOneDThreeDTestCCProblem, INHERITS_FROM(CCModel, OnePOneDThreeDTestProblem));
SET_TYPE_PROP(OnePOneDThreeDTestProblem, Grid, Dune::FoamGrid<1, 3>);

NEW_TYPE_TAG(OnePTwoDThreeDTestProblem, INHERITS_FROM(OnePDupuitThiemProblem));
NEW_TYPE_TAG(OnePTwoDThreeDTestBoxProblem, INHERITS_FROM(BoxModel, OnePTwoDThreeDTestProblem));
NEW_TYPE_TAG(OnePTwoDThreeDTestCCProblem, INHERITS_FROM(CCModel, OnePTwoDThreeDTestProblem));
SET_TYPE_PROP(OnePTwoDThreeDTestProblem, Grid, Dune::FoamGrid<2, 3>);
#endif

// Enable gravity
SET_BOOL_PROP(OnePDupuitThiemProblem, ProblemEnableGravity, false);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * This problem uses an extrusion factor which allows the 2D - simulation of a complete radial domain whith an
 * extraction well in its center.
 * Modelled is just a lateral section of the domain which is than extruded to 3D using above-mentioned extrusion factor.
 *
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the left boundary (Neumann), where water is extracted by a well
 * and the right boundary, where a Dirichlet BC for pressure holds.
 * .
 *
 * The domain has a constant permeability (\f$K=10e-10\f$).
 */
template <class TypeTag>
class OnePDupuitThiemProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox)
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;


    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePDupuitThiemProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView), eps_(1e-6)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        qWell_ = 1e-3; // well pumping rate in m^3/s
    }

    void postTimeStep()
    {
        const SolutionVector& solution = this->model().curSol();

        //calculate the discrete L2-Norm
        Scalar LTwoNorm = 0.0;

        //get the Gaussian quadrature rule for intervals

        for (const auto& element : elements(this->gridView()))
        {
            const unsigned int eIdx = this->elementMapper().index(element);
            const auto geometry = element.geometry();
            const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(geometry.type(), 1/*order*/);
            for(auto&& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const Scalar pe = pressureExact(globalPos);
                const Scalar integrationElement = geometry.integrationElement(qp.position());
                Scalar p = 0.0;
                if (!isBox)
                    p = solution[eIdx][pressureIdx];
                else
                {
                    // do interpolation with ansatz functions
//                     std::vector<Dune::FieldVector<Scalar, 1> > shapeValues;
//                     const auto& localFiniteElement = feCache_.get(geometry.type());
//                     localFiniteElement.localBasis().evaluateFunction(qp.position(), shapeValues);
//                     for (unsigned int i = 0; i < shapeValues.size(); ++i)
//                         p += shapeValues[i]*solution[this->model().dofMapper().subIndex(element, i, dim)][pressureIdx];
                }
                LTwoNorm += (p - pe)*(p - pe)*qp.weight()*integrationElement;
            }
        }
        LTwoNorm = std::sqrt(LTwoNorm);

        logFile_.open(this->name() + ".log", std::ios::app);
        logFile_ << "[ConvergenceTest] L2-norm(pressure) = " << LTwoNorm  << " hMax = " << this->bBoxMax()[0]/this->gridView().size(0) << std::endl;
        logFile_.close();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
        {
            //Here we calculate the analytical solution
            typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
            const unsigned numDofs = this->model().numDofs();

           //create required scalar fields
           ScalarField *pAnalytical = this->resultWriter().allocateManagedBuffer(numDofs);

           for (const auto& element : elements(this->gridView()))
            {
                const auto& fvGeometry = this->model().fvGeometries(element);
                for (auto&& scv : fvGeometry.scvs() )
                {
                    const auto globalPos = scv.dofPosition();
                    const auto globalIdx = scv.dofIndex();
                    (*pAnalytical)[globalIdx] = pressureExact(globalPos);
                }
            }
            this->resultWriter().attachDofData(*pAnalytical, "pA", isBox);
        }

    /*!
     * \brief Return the pressure at a given position computed with the Dupuit-Thiem
     *        analytical solution for confined aquifers under steady state conditions
     */
    Scalar pressureExact(const GlobalPosition& globalPos) const
    {
        const Scalar pRight = dirichletAtPos(this->bBoxMax())[pressureIdx];
        const Scalar k = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);
        const Scalar m = this->bBoxMax()[1]-this->bBoxMin()[1];
        const Scalar densityW = 1e3;
        const Scalar mobilityW = 1/1e-3;
        const Scalar h0 = pRight/(9.81*densityW);
        const Scalar conductivity = k * mobilityW * densityW * 9.81;
        const Scalar drawdown = qWell_ / (2*M_PI*m*conductivity)*std::log((this->bBoxMax()[0])/(globalPos[0]));
        const Scalar h = h0-drawdown;
        const Scalar p = h*9.81*1000;
        return p;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     *//*PrimaryVariables*/
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0);
    }
    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        Scalar eps = 1.0e-6;
        if (globalPos[0] > this->bBoxMax()[0]-eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0);
        values[pressureIdx] = 1.0e+5;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    PrimaryVariables neumann(const Element &element, const SubControlVolumeFace &scvFace) const
    {
        PrimaryVariables priVars(0);
        auto globalPos = scvFace.center();
        if(globalPos[0] < this->bBoxMin()[0] + 1e-6)
            priVars[conti0EqIdx] = qWell_ / this->bBoxMax()[1] / extrusionFactorAtPos(globalPos) * 1000 ;
        return priVars;
    }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     *
     * In this problem, the extrusion factor equals the circumference of a cylinder at given radius
     *
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        return  2*M_PI*(globalPos[0]);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const SubControlVolume& scv) const
    {
        PrimaryVariables priVars(0);
        priVars[pressureIdx] = 1.0e+5;
        return priVars;
    }

    // \}

private:
    std::string name_;
    Scalar eps_;
    Scalar qWell_;
    std::ofstream logFile_;
};
} //end namespace

#endif
