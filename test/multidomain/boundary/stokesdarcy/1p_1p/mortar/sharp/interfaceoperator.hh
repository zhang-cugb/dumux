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
#ifndef DUMUX_MORTAR_ONEP_SHARP_INTERFACE_OPERATOR_HH
#define DUMUX_MORTAR_ONEP_SHARP_INTERFACE_OPERATOR_HH

#include <dune/common/indices.hh>
#include <dune/istl/operators.hh>

#include <dumux/discretization/projection/projector.hh>
#include <dumux/discretization/functionspacebasis.hh>
#include <dumux/multidomain/glue.hh>

#include "../reconstructionhelper.hh"
#include "projector.hh"

namespace Dumux {

//! HACK for staggered. Provide FunctionSpaceBasisTraits.
//! We use piecewise constants here.
template< class GridGeometry >
struct FunctionSpaceBasisTraits<GridGeometry, DiscretizationMethod::staggered>
{ using GlobalBasis = Dune::Functions::LagrangeBasis<typename GridGeometry::GridView, /*order*/0>; };

template<class SolutionVector,
         class MortarGridGeometry,
         class Solver1, class Reconstructor1,
         class Solver2, class Reconstructor2>
class SharpMortarInterfaceOperator
: public Dune::LinearOperator< SolutionVector, SolutionVector >
{
    using FieldType =  typename SolutionVector::field_type;
    using Projector = SharpMortarProjector<SolutionVector>;

    using GridGeometry1 = typename Solver1::FVGridGeometry;
    using GridGeometry2 = typename Solver2::FVGridGeometry;
    using CoupledScvfMap1 = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry1>;
    using CoupledScvfMap2 = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry2>;

public:
    explicit SharpMortarInterfaceOperator(std::shared_ptr<Solver1> solver1,
                                          std::shared_ptr<Solver2> solver2,
                                          std::shared_ptr<MortarGridGeometry> mortarGG)
    : solver1_(solver1)
    , solver2_(solver2)
    , mortarGridGeometry_(mortarGG)
    {
        using Helper = MortarReconstructionHelper;
        coupledScvfMap1_ = Helper::findCoupledScvfs(*solver1->gridGeometryPointer(), *mortarGG);
        coupledScvfMap2_ = Helper::findCoupledScvfs(*solver2->gridGeometryPointer(), *mortarGG);

        const auto& gg1 = *solver1->gridGeometryPointer();
        const auto& gg2 = *solver2->gridGeometryPointer();

        const auto glue1 = makeGlue(gg1, *mortarGG);
        const auto glue2 = makeGlue(gg2, *mortarGG);

        const auto projector1 = makeProjectorPair(getFunctionSpaceBasis(gg1), getFunctionSpaceBasis(*mortarGG), glue1).second;
        const auto projector2 = makeProjectorPair(getFunctionSpaceBasis(gg2), getFunctionSpaceBasis(*mortarGG), glue2).second;

        projector_ = std::make_shared<Projector>(projector1, projector2);
    }

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply(const SolutionVector& x, SolutionVector& r) const
    {
        using namespace Dune::Indices;

        // project mortar pressure into sub-domains
        auto x1 = x;
        auto x2 = x;

        const auto isNegative1 = solver1_->problemPointer()->isOnNegativeMortarSide();
        const auto isNegative2 = solver2_->problemPointer()->isOnNegativeMortarSide();

        if (isNegative1) x1 *= -1.0;
        if (isNegative2) x2 *= -1.0;

        auto flux = projector_->projectMortarToSubDomain(x1, x2);

        // sub-problem negate the fluxes already (avoid doing it twice)
        if (isNegative1) flux[_0] *= -1.0;
        if (isNegative2) flux[_1] *= -1.0;

        solver1_->problemPointer()->setMortarProjection(flux[_0]);
        solver2_->problemPointer()->setMortarProjection(flux[_1]);

        solver1_->solve();
        solver2_->solve();

        using R1 = Reconstructor1;
        using R2 = Reconstructor2;

        auto pressure1 = R1::template recoverSolution<SolutionVector>(*solver1_->gridGeometryPointer(),
                                                                            *solver1_->gridVariablesPointer(),
                                                                            *solver1_->solutionPointer(),
                                                                            coupledScvfMap1_);

        auto pressure2 = R2::template recoverSolution<SolutionVector>(*solver2_->gridGeometryPointer(),
                                                                            *solver2_->gridVariablesPointer(),
                                                                            *solver2_->solutionPointer(),
                                                                            coupledScvfMap2_);

        auto pressureTuple = projector_->projectSubDomainToMortar(pressure1, pressure2);
        if (isNegative1) pressureTuple[_0] *= -1.0;
        if (isNegative2) pressureTuple[_1] *= -1.0;

        r = pressureTuple[_0];
        r += pressureTuple[_1];
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd(FieldType alpha, const SolutionVector& x, SolutionVector& y) const
    {
        SolutionVector yTmp;

        apply(x, yTmp);
        yTmp *= alpha;

        y += yTmp;
    }

    //! Category of the solver (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

private:
    std::shared_ptr<Solver1> solver1_;
    std::shared_ptr<Solver2> solver2_;
    std::shared_ptr<MortarGridGeometry> mortarGridGeometry_;

    std::shared_ptr<Projector> projector_;

    CoupledScvfMap1 coupledScvfMap1_;
    CoupledScvfMap2 coupledScvfMap2_;
};

} // end namespace Dumux

#endif
