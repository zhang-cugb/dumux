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
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_STOKES_DARCY_MORTAR_SHARP_PRECONDITIONER_HH
#define DUMUX_STOKES_DARCY_MORTAR_SHARP_PRECONDITIONER_HH

#include <dune/istl/preconditioner.hh>
#include "interfaceoperator.hh"

namespace Dumux {

template<class SolutionVector,
         class MortarGridGeometry,
         class Solver1, class Reconstructor1,
         class Solver2, class Reconstructor2 >
class SharpMortarPreconditioner
: public Dune::Preconditioner<SolutionVector, SolutionVector>
{
    using FieldType =  typename SolutionVector::field_type;

    using GridGeometry1 = typename Solver1::FVGridGeometry;
    using GridGeometry2 = typename Solver2::FVGridGeometry;
    using CoupledScvfMap1 = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry1>;
    using CoupledScvfMap2 = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry2>;

public:
    explicit SharpMortarPreconditioner(std::shared_ptr<Solver1> solver1,
                                      std::shared_ptr<Solver2> solver2,
                                      std::shared_ptr<MortarGridGeometry> mortarGG)
    : solver1_(solver1)
    , solver2_(solver2)
    , mortarGridGeometry_(mortarGG)
    {}

    /*!
     * \brief Prepare the preconditioner.
     */
    virtual void pre (SolutionVector& x, SolutionVector& b) {}

    /*!
     * \brief Apply one step of the preconditioner to the system A(v)=d.
     */
    virtual void apply (SolutionVector& r, const SolutionVector& x)
    {
        r = x;
    }

    /*!
     * \brief Clean up.
     */
    virtual void post (SolutionVector& x) {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

private:
    std::shared_ptr<Solver1> solver1_;
    std::shared_ptr<Solver2> solver2_;
    std::shared_ptr<MortarGridGeometry> mortarGridGeometry_;
};

} // end namespace Dumux

#endif
