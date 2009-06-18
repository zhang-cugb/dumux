// $Id:$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/**
 * @file
 * @brief  Definition of a problem, where air is injected under a low permeable layer
 * @author Bernd Flemisch, Klaus Mosthaf
 */
#ifndef DUNE_INJECTIONPROBLEM_HH
#define DUNE_INJECTIONPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluids/water_air.hh>

#include <dumux/material/multicomponentrelations.hh>

#include <dumux/boxmodels/2p2c/2p2cboxmodel.hh>

#include "injectionsoil.hh"


namespace Dune
{

template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem, INHERITS_FROM(BoxTwoPTwoC));

// Set the grid type
SET_PROP(InjectionProblem, Grid)
{
#if ENABLE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::SGrid<2, 2> type;
    //typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(InjectionProblem, Problem)
{
    typedef Dune::InjectionProblem<TTAG(InjectionProblem)> type;
};

// Set the wetting phase
SET_TYPE_PROP(InjectionProblem, WettingPhase, Dune::Liq_WaterAir);

// Set the non-wetting phase
SET_TYPE_PROP(InjectionProblem, NonwettingPhase, Dune::Gas_WaterAir);

// Set multi-component relations
SET_TYPE_PROP(InjectionProblem, MultiComp, Dune::CWaterAir)

// Set the soil properties
SET_PROP(InjectionProblem, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    
public:
    typedef Dune::InjectionSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem, EnableGravity, true);
}


/*! 
 * \brief Problem where air is injected under a low permeable layer.
 *
 * Air enters the domain at the right boundary and migrates upwards.
 * Problem was set up using the rect2d.dgf grid.
 */
template <class TypeTag = TTAG(InjectionProblem) >
class InjectionProblem : public TwoPTwoCBoxProblem<TypeTag, InjectionProblem<TypeTag> >
{
    typedef InjectionProblem<TypeTag>             ThisType;
    typedef TwoPTwoCBoxProblem<TypeTag, ThisType> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        // Grid and world dimension
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;

    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef typename GridView::template Codim<dim>::Entity      Vertex;
    typedef typename GridView::IntersectionIterator             IntersectionIterator;
  
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

public:
    InjectionProblem(const GridView &gridView)
        : ParentType(gridView)
    {
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
    const char *name() const
    { return "injection"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 283.15; // -> 10°C
    };
    
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*! 
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        if (globalPos[0] < eps_)
            values = BoundaryConditions::dirichlet;
        else
            values = BoundaryConditions::neumann;
    }

    /*! 
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVarVector           &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

    /*! 
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVarVector           &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        values = 0;
        if (globalPos[1] < 15 && globalPos[1] > 5) {
            values[Indices::comp2Mass(Indices::nComp)] = -1e-3;
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*! 
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*! 
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVarVector        &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

    /*! 
     * \brief Return the initial phase state inside a control volume.
     */
    int initialPhaseState(const Vertex       &vert,
                          int                &globalIdx,
                          const GlobalPosition &globalPos) const
    { return Indices::wPhaseOnly; }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVarVector       &values,
                  const GlobalPosition   &globalPos) const
    {
        Scalar densityW = this->wettingPhase().density(temperature(), 1e5);

        values[Indices::pW] = 1e5 - densityW*this->gravity()[1]*(depthBOR_ - globalPos[1]);
        values[Indices::sNorX] = 0;
    }

    static const Scalar depthBOR_ = 800.0; // [m]
    static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
