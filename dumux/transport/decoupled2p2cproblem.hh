// $Id:$

/*****************************************************************************
* Copyright (C) 2008 by Jochen Fritz                                         *
* Institute of Hydraulic Engineering                                         *
* University of Stuttgart, Germany                                           *
* email: <givenname>.<name>@iws.uni-stuttgart.de                             *
*                                                                            *
* This program is free software; you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation; either version 2 of the License, or          *
* (at your option) any later version, as long as this copyright notice       *
* is included in its original form.                                          *
*                                                                            *
* This program is distributed WITHOUT ANY WARRANTY.                          *
*****************************************************************************/

#ifndef DECOUPLED2P2CPROBLEM_HH
#define DECOUPLED2P2CPROBLEM_HH

#include <iostream>
#include <iomanip>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/disc/operators/boundaryconditions.hh>

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/operators/boundaryconditions2p2c.hh>
#include <dumux/fractionalflow/variableclass2p2c.hh>



namespace Dune
{
/**
 * \ingroup decoupled2p2cproblems
 *
 * \brief Base class for the definition of 2p2c problems
 *
 * This base class defines all boundary and initial functions which are needed
 * for a decoupled 2p2c computation.
 */
template<class GridView, class Scalar>
class DecoupledProblem2p2c
{
    enum {dim=GridView::dimension};
    typedef typename GridView::Traits::template Codim<0>::Entity Entity;

public:
    //! Type of transport boundary condition.
    /**    either the concentration or the saturation have to be defined
     * on boundaries with dirichlet pressure BCs.
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual BoundaryConditions2p2c::Flags bc_type (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                                   const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Type of transport initial condition.
    /**    either the concentration or the saturation have to be defined
     * as initial condition.
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual BoundaryConditions2p2c::Flags initcond_type (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                                         const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Type of pressure boundary condition.
    /**    Pressure (dirichlet) or flux (neumann) have to be defined on boundaries.
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual BoundaryConditions::Flags press_bc_type (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                                     const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Feed mass fraction boundary condition (dimensionless)
    /** Feed concentration is the (global) mass fraction (dimensionless) of component 1 in the mixture
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar dirichletConcentration (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                           const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Saturation boundary condition (dimensionless)
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar dirichletSat (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                 const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Pressure (dirichlet) boundary condition \f$ [Pa] \f$
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar dirichlet (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                              const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Flux (neumann) boundary condition \f$ [\frac{kg}{m^2 \cdot s}] \f$
    /** Returns a vector with one entry for each component. The flux is expressed in \f$ [\frac{kg}{m^2 \cdot s}] \f$
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual FieldVector<Scalar,2> neumann (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                           const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Source of components \f$ [\frac{kg}{m^3 \cdot s}] \f$
    /** Returns a vector with one entry for each component.
     * Describes the source of the components per unit area in \f$ [\frac{kg}{m^3 \cdot s}] \f$
     * @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual FieldVector<Scalar,2> source (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                          const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Saturation initial condition (dimensionless)
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar initSat (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                            const FieldVector<Scalar,dim>& localPos) const = 0;

    //! Feed mass fraction initial condition (dimensionless)
    /** @param globalPos global coordinates
     * @param element reference to the cell for which the function is to be evaluated
     * @param localPos local coordinates inside element
     */
    virtual Scalar initConcentration (const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                      const FieldVector<Scalar,dim>& localPos) const = 0;

    //! gravity vector \f$ [\frac{m}{s^2}] \f$
    virtual const FieldVector<Scalar,dim> gravity()
    {
        FieldVector<Scalar,dim> gravity_(0.);
        return gravity_;
    }

    //! Temperature \f$ [K] \f$
    virtual Scalar temp(const FieldVector<Scalar,dim>& globalPos, const Entity& element,
                                      const FieldVector<Scalar,dim>& localPos) = 0;

    //! Constructor
    /**
     * \param var a VariableClass2p2c object
     * \param liq class defining the properties of the liquid phase
     * \param gas class defining the properties of the gaseous phase
     * \param s the soil properties
     * \param law the relative permeability and capillary pressure definintions
     */
    DecoupledProblem2p2c(Dune::VariableClass2p2c<GridView, Scalar>& var, Liquid_GL& liq, Gas_GL& gas, Matrix2p<typename GridView::Grid, Scalar>& s,
                         TwoPhaseRelations<typename GridView::Grid, Scalar>& law = *(new TwoPhaseRelations<typename GridView::Grid,Scalar>),const bool cap = false)
        :variables(var), liquidPhase(liq), gasPhase(gas), soil(s), materialLaw(law), capillary(cap)
    {
    }

    virtual ~DecoupledProblem2p2c ()
    {
    }

    VariableClass2p2c<GridView, Scalar>& variables;
    Liquid_GL& liquidPhase;
    Gas_GL& gasPhase;
    Matrix2p<typename GridView::Grid, Scalar>& soil;
    TwoPhaseRelations<typename GridView::Grid, Scalar>& materialLaw;
    const bool capillary;
};

}
#endif
