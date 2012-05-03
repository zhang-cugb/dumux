// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *                 2012 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for all problems which use the two-phase box model
 */
#ifndef DUMUX_POROUS_MEDIA_BOX_PROBLEM_HH
#define DUMUX_POROUS_MEDIA_BOX_PROBLEM_HH

#include "boxproperties.hh"

#include <dumux/boxmodels/common/boxproblem.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters object
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
}

/*!
 * \ingroup BoxBaseProblems
 * \brief Base class for all porous media box problems
 */
template<class TypeTag>
class PorousMediaBoxProblem : public BoxProblem<TypeTag>
{
    typedef BoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param verbose Turn verbosity on or off
     */
    PorousMediaBoxProblem(TimeManager &timeManager,
                const GridView &gridView,
                const bool verbose = true)
        : ParentType(timeManager, gridView),
          gravity_(0)
    {
        newSpatialParams_ = true;
        spatialParams_ = new SpatialParams(gridView);

        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim-1]  = -9.81;
    }

    ~PorousMediaBoxProblem()
    {
        if (newSpatialParams_)
            delete spatialParams_;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * This is the discretization specific interface for the box
     * method. By default it just calls temperature(pos).
     *
     * \param element The DUNE Codim<0> enitiy which intersects with
     *                the finite volume.
     * \param fvGeometry The finite volume geometry of the element.
     * \param scvIdx The local index of the sub control volume inside the element
     */
    Scalar boxTemperature(const Element &element,
                          const FVElementGeometry fvGeometry,
                          const int scvIdx) const
    { return asImp_().temperatureAtPos(fvGeometry.subContVol[scvIdx].global); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param pos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &pos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); };

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is the box discretization specific interface. By default
     * it just calls gravityAtPos().
     */
    const DimVector &boxGravity(const Element &element,
                                const FVElementGeometry &fvGeometry,
                                const int scvIdx) const
    { return asImp_().gravityAtPos(fvGeometry.subContVol[scvIdx].global); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is discretization independent interface. By default it
     * just calls gravity().
     */
    const DimVector &gravityAtPos(const GlobalPosition &pos) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>EnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const DimVector &gravity() const
    { return gravity_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return *spatialParams_; }
    
    DUMUX_DEPRECATED_MSG("use spatialParams instead")
    SpatialParams &spatialParameters()
    { return spatialParams(); }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

    DUMUX_DEPRECATED_MSG("use spatialParams instead")
    const SpatialParams &spatialParameters() const
    { return spatialParams(); }

    // \}

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    DimVector gravity_;

    // fluids and material properties
    SpatialParams*  spatialParams_;
    bool newSpatialParams_;
};

}

#endif
