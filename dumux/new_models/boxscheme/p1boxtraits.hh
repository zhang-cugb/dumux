/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_PWSN_BOX_TRAITS_HH
#define DUMUX_PWSN_BOX_TRAITS_HH

#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dumux/fvgeometry/fvelementgeometry.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxScheme.
 */
namespace Dune
{
    /*!
     * \brief Specify the shape functions, operator assemblers, etc
     *        used for the BoxScheme.
     */
    template<class DT, class Grid, int numEqns>
    class P1BoxTraits
    {
    private:
        typedef typename Grid::ctype     CoordScalar;

        enum {
            GridDim = Grid::dimension,
            WorldDim = Grid::dimensionworld,
        };
    
    public:
        enum {
            numEq = numEqns
        };
        
        //! A vector of all primary variables at a point
        typedef Dune::FieldVector<DT, numEq> UnknownsVector;
        //! boundary condition vector
        typedef Dune::FieldVector<Dune::BoundaryConditions::Flags, 
                                  numEq>         BoundaryTypeVector;
        
        /*!
         * \brief Represents a local spatial function.
         *
         * A field vector is attached at each vertex of the cell.
         */
        typedef BlockVector<UnknownsVector> LocalFunction;
        
        
        //! The finite volume cell segments within a finite element cell
        typedef Dune::FVElementGeometry<Grid>  FVElementGeometry;
        
        //! a single of shape function used for the BoxFunction inside
        //! cells.
        typedef Dune::LagrangeShapeFunctionSetContainer<CoordScalar,
                                                        DT,
                                                        GridDim> ShapeFunctionSetContainer;
        
        //! The actual shape functions which are being used. If a 
        //! grid only contains simplices or tetrahedra, it is more
        //! efficent to use LagrangeShapeFunctions::p1cube, or
        //! LagrangeShapeFunctions::p1simplex
        //!
        //! TODO: Use specialization to take advantage of simplex grids
        //!       and structured grids.
        static const ShapeFunctionSetContainer &shapeFunctions;
        
        //! The function which represents a solution for a fixed time
        //! step. We use first-order vertex centered FE polynomials.
        typedef Dune::LeafP1Function<Grid, 
                                     DT, 
                                     numEq>   SpatialFunction;
        
        //! The OperatorAssembler which assembles the global stiffness
        //! matrix
        typedef Dune::LeafP1OperatorAssembler<Grid,
                                              DT,
                                              numEq>  JacobianAssembler;
    };

    // this is butt-ugly, but the only way I could come up with to
    // initialize a static const member of a class
    template<class DT, class Grid, int PrimaryVars>
    const typename P1BoxTraits<DT, Grid, PrimaryVars>::ShapeFunctionSetContainer &
      P1BoxTraits<DT, Grid, PrimaryVars>::shapeFunctions
         =  Dune::LagrangeShapeFunctions<typename Grid::ctype, 
                                         DT,
                                         Grid::dimension>::general;

}

#endif
