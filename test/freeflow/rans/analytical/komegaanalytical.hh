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
*
* \brief A analytic 1/2-D solution for thekomega turbulence model.
*
* The 1-D analytic solution is given by
* \f[ v_\text{x} = 2 \cdot x^3 \f]
* \f[ p = 2 - 2 \cdot x \f]
* \f[ k = 1 + x^2 \f]
* \f[ \omega = 2 - x^2\f]
*
* The 2-D analytic solution is given by
* \f[ v_\text{x} = 2 \cdot x^3 \f]
* \f[ v_\text{y} = 2 - 2 \cdot y^3 \f]
* \f[ p = 2 - 2 \cdot x y \f]
* \f[ k = 1 + x^2 y^2 \f]
* \f[ \omega = 2 - x^2 y^2 \f]
*/

#ifndef DUMUX_FREEFLOW_RANS_KOMEGA_ANALYTICAL_HH
#define DUMUX_FREEFLOW_RANS_KOMEGA_ANALYTICAL_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {

template<class Scalar, class GlobalPosition, class PrimaryVariables, class ModelTraits, class TurbulentConstants, int dim>
class KOmegaAnalytical
{
public:
    using VelocityVector = Dune::FieldVector<Scalar, dim>;
    using VelocityGradientMatrix = Dune::FieldVector<VelocityVector, dim>;
    using ScalarGradientVector = Dune::FieldVector<Scalar, dim>;

    using Indices = typename ModelTraits::Indices;
    using GravityVector = Dune::FieldVector<Scalar, dim>;

    KOmegaAnalytical(const Scalar density, const Scalar kinematicViscosity, const TurbulentConstants turbConstants)
    : density_(density)
    , kinematicViscosity_(kinematicViscosity)
    , turbConstants_(turbConstants)
    {}

    /*!
    * \brief Returns the analytical solution of the problem at a given time and position.
    *
    * \param globalPos The global position
    * \param time The current simulation time
    */
    PrimaryVariables solution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        values[Indices::pressureIdx] = sourceMassBalanceAtPos(globalPos);

        const VelocityVector momentumSolution = sourceMomentumBalanceAtPos(globalPos);
        values[Indices::momentumXBalanceIdx] = momentumSolution[0];
        if (dim > 1)
            values[Indices::momentumYBalanceIdx] = momentumSolution[1];

        values[Indices::turbulentKineticEnergyIdx] = sourceTurbulentKineticEnergyBalanceAtPos(globalPos);
        values[Indices::dissipationIdx] = sourceDissipationBalanceAtPos(globalPos);

        return values;
    }

private:

    Scalar pressure(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            return 2.0 - 2.0 * x * y;
        }
        else
            return 2.0 - 2.0 * x;
    }

    VelocityVector velocity(const GlobalPosition& globalPos) const
    {
        VelocityVector v(0.0);
        const Scalar x = globalPos[0];
        v[0] = 2.0 * x * x * x;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            v[1] = 2.0 - 2.0 * y * y * y;
        }
        return v;
    }

    Scalar turbulentKineticEnergy(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            return 1.0 + x * x * y * y;
        }
        else
            return 1.0 + x * x;
    }

    Scalar dissipation(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            return 2.0 - x * x * y * y;
        }
        else
            return 2.0 - x * x;
    }

    Scalar limitedDissipation(const GlobalPosition& globalPos) const
    {
        // use the Dissipation limiter proposed in wilcox2008
        using std::sqrt;
        using std::max;
        return (7.0 / 8.0) * sqrt(2.0 * shearStressTensorProduct(globalPos) / turbConstants_.betaK);
    }

    //! \brief The gradients
    VelocityGradientMatrix dvdx(const GlobalPosition& globalPos) const
    {
        VelocityGradientMatrix dvdx;
        const Scalar x = globalPos[0];
        dvdx[0][0] = 6.0 * x * x;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dvdx[0][1] = 0.0;
            dvdx[1][0] = 0.0;
            dvdx[1][1] = -6.0 * y * y;
        }
        return dvdx;
    }

    VelocityGradientMatrix dv2dx(const GlobalPosition& globalPos) const
    {
        VelocityGradientMatrix dv2dx;
        for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
        {
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                dv2dx[velIdx][dimIdx] = dvdx(globalPos)[velIdx][dimIdx] * velocity(globalPos)[dimIdx]
                                      + dvdx(globalPos)[dimIdx][dimIdx] * velocity(globalPos)[velIdx];
            }
        }

        return dv2dx;
    }

    VelocityGradientMatrix dvdx2(const GlobalPosition& globalPos) const
    {
        VelocityGradientMatrix dvdx2;
        const Scalar x = globalPos[0];
        dvdx2[0][0] = 12.0 * x;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dvdx2[0][1] = 0.0;
            dvdx2[1][0] = 0.0;
            dvdx2[1][1] = -12.0 * y;
        }
        return dvdx2;
    }

    ScalarGradientVector dpdx(const GlobalPosition& globalPos) const
    {
        ScalarGradientVector dpdx(0.0);
        const Scalar x = globalPos[0];
        dpdx[0] = -2.0;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dpdx[0] = -2.0 * y;
            dpdx[1] = -2.0 * x;
        }
        return dpdx;
    }

    ScalarGradientVector dkdx(const GlobalPosition& globalPos) const
    {
        ScalarGradientVector dkdx(0.0);
        const Scalar x = globalPos[0];
        dkdx[0] = 2.0 * x;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dkdx[0] = 2.0 * x * y * y;
            dkdx[1] = 2.0 * y * x * x;
        }
        return dkdx;
    }

    ScalarGradientVector dkdx2(const GlobalPosition& globalPos) const
    {
        ScalarGradientVector dkdx2(0.0);
        const Scalar x = globalPos[0];
        dkdx2[0] = 2.0;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dkdx2[0] = 2.0 * y * y;
            dkdx2[1] = 2.0 * x * x;
        }
        return dkdx2;
    }

    ScalarGradientVector dwdx(const GlobalPosition& globalPos) const
    {
        ScalarGradientVector dwdx(0.0);
        const Scalar x = globalPos[0];
        dwdx[0] = -2.0 * x;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dwdx[0] = -2.0 * x * y * y;
            dwdx[1] = -2.0 * x * x * y;
        }
        return dwdx;
    }

    ScalarGradientVector dwdx2(const GlobalPosition& globalPos) const
    {
        ScalarGradientVector dwdx2(0.0);
        const Scalar x = globalPos[0];
        dwdx2[0] = -2.0;
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dwdx2[0] = -2.0 * y * y;
            dwdx2[1] = -2.0 * x * x;
        }
        return dwdx2;
    }

    // eddy viscosity and its gradient
    Scalar nut(const GlobalPosition& globalPos) const
    {
        using std::max;
        // use the Dissipation limiter proposed in wilcox2008
        static const auto enableKOmegaDissipationLimiter
            = getParam<bool>("Problem.KOmega.EnableDissipationLimiter", true);
        if (enableKOmegaDissipationLimiter)
            return turbulentKineticEnergy(globalPos) / max(dissipation(globalPos), limitedDissipation(globalPos));

        return turbulentKineticEnergy(globalPos) / dissipation(globalPos);
    }

    ScalarGradientVector dnutdx(const GlobalPosition& globalPos) const
    {
        using std::sqrt;
        using std::pow;
        ScalarGradientVector dnutdx(0.0);
        const Scalar x = globalPos[0];

        static const auto enableKOmegaDissipationLimiter
            = getParam<bool>("Problem.KOmega.EnableDissipationLimiter", true);
        if (enableKOmegaDissipationLimiter)
        {    const bool inLimiterSpace = limitedDissipation(globalPos) > dissipation(globalPos);
            if (inLimiterSpace)
            {
                dnutdx[0] = (4. * sqrt(2.)) / ( 21. * x * sqrt(x*x*x*x / turbConstants_.betaK));
                if (dim > 1)
                {
                    const Scalar y = globalPos[1];
                    dnutdx[0] = (4. * sqrt(2.) * x * (x*x - y*y*y*y*y*y))
                              / ( 21. * turbConstants_.betaK * pow((x*x*x*x + y*y*y*y) / turbConstants_.betaK, (3./2.)));

                    dnutdx[1] = (4. * sqrt(2.) * y * (x*x*x*x*x*x - y*y))
                              / ( 21. * turbConstants_.betaK * pow((x*x*x*x + y*y*y*y) / turbConstants_.betaK, (3./2.)));
                }
            }
        }

        dnutdx[0] = (6*x) / ((2-(x*x)) * (2-(x*x)));
        if (dim > 1)
        {
            const Scalar y = globalPos[1];
            dnutdx[0] = 6*x*y*y
                    / ((2-(x*x*y*y)) * (2-(x*x*y*y)));
            dnutdx[1] = 6*x*x*y
                    / ((2-(x*x*y*y)) * (2-(x*x*y*y)));
        }
        return dnutdx;
    }

    // shearStressTensorProduct
    Scalar shearStressTensorProduct(const GlobalPosition& globalPos) const
    {
        Scalar shearStressTensorProduct = 0.0;
        for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
        {
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                shearStressTensorProduct += 0.5 * (dvdx(globalPos)[velIdx][dimIdx] + dvdx(globalPos)[dimIdx][velIdx])
                                          * 0.5 * (dvdx(globalPos)[velIdx][dimIdx] + dvdx(globalPos)[dimIdx][velIdx]);
            }
        }
        // Should be 36x^4 +36y^4 (1D: 36x^4)
        return shearStressTensorProduct;
    }

    //! \brief Source term values
    Scalar sourceMassBalanceAtPos(const GlobalPosition& globalPos) const
    {
        // term div(rho*v)
        Scalar divV(0.0);
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            divV += dvdx(globalPos)[dimIdx][dimIdx];
        }
        return density_ * divV;
    }

    VelocityVector sourceMomentumBalanceAtPos(const GlobalPosition& globalPos) const
    {
        VelocityVector momentum(0.0);
        GravityVector gravity(0.0);
        const bool enableGravity = getParam<bool>("Problem.EnableGravity");
        if (enableGravity)
            gravity[dim-1] = -9.81;

        for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
        {
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                // inertia term
                momentum[velIdx] += density_ * dv2dx(globalPos)[velIdx][dimIdx];

                // viscous term (molecular)
                momentum[velIdx] -= density_ * kinematicViscosity_ * dvdx2(globalPos)[velIdx][dimIdx];
                momentum[velIdx] -= density_ * kinematicViscosity_ * dvdx2(globalPos)[dimIdx][velIdx];

                // viscous term (turbulent) (using product rule)
                momentum[velIdx] -= density_ * ((dnutdx(globalPos)[velIdx] * dvdx(globalPos)[velIdx][dimIdx])
                                              + (nut(globalPos) * dvdx2(globalPos)[velIdx][dimIdx]));
                momentum[velIdx] -= density_ * ((dnutdx(globalPos)[dimIdx] * dvdx(globalPos)[velIdx][dimIdx])
                                              + (nut(globalPos) * dvdx2(globalPos)[dimIdx][velIdx]));
            }

            // pressure term
            momentum[velIdx] += dpdx(globalPos)[velIdx];

            // gravity term
            if (enableGravity)
                momentum[velIdx] += density_ * gravity[velIdx];
        }
        return momentum;
    }

    Scalar sourceTurbulentKineticEnergyBalanceAtPos(const GlobalPosition& globalPos) const
    {
        Scalar turbulentKineticEnergySol = 0.0;

        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            // advection term (using product rule)
            turbulentKineticEnergySol += density_ * ( (dvdx(globalPos)[dimIdx][dimIdx] * turbulentKineticEnergy(globalPos))
                                                    + (velocity(globalPos)[dimIdx] * dkdx(globalPos)[dimIdx]) );
            // Molecular diffusion term
            turbulentKineticEnergySol -= kinematicViscosity_ * dkdx2(globalPos)[dimIdx];

            // Turbulent diffusion term (using product rule)
            turbulentKineticEnergySol -= density_ *  turbConstants_.sigmaK * ( (dnutdx(globalPos)[dimIdx] * dkdx(globalPos)[dimIdx])
                                                                             + (nut(globalPos) * dkdx2(globalPos)[dimIdx]));
        }

        // Production term
        static const auto enableKOmegaProductionLimiter = getParam<bool>("Problem.KOmega.EnableProductionLimiter", false);
        const Scalar kProductionOriginal= 2.0 * nut(globalPos) * shearStressTensorProduct(globalPos);
        Scalar kProductionAlternative = 0.0;
        if (enableKOmegaProductionLimiter)
        {
            // Production limiter (2006 K-limiter via P term)
            kProductionAlternative += 20.0 * turbConstants_.betaK * turbulentKineticEnergy(globalPos) * dissipation(globalPos);
        }
        turbulentKineticEnergySol -= std::min(kProductionOriginal, kProductionAlternative);

        // Destruction term
        turbulentKineticEnergySol += turbConstants_.betaK * density_ * turbulentKineticEnergy(globalPos) * dissipation(globalPos);

        return turbulentKineticEnergySol;
    }

    Scalar sourceDissipationBalanceAtPos(const GlobalPosition& globalPos) const
    {
        Scalar dissipationSol = 0.0;
        // Advection & Diffusion (using product rule)
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            // Advection term
            dissipationSol += ( (dvdx(globalPos)[dimIdx][dimIdx] * dissipation(globalPos))
                              + (velocity(globalPos)[dimIdx] * dwdx(globalPos)[dimIdx]) );

            // Molecular diffusion term
            dissipationSol -= kinematicViscosity_ * dwdx2(globalPos)[dimIdx];

            // Turbulent diffusion term
            dissipationSol -= turbConstants_.sigmaOmega * ( (dnutdx(globalPos)[dimIdx] * dwdx(globalPos)[dimIdx])
                                                           + (nut(globalPos) * dwdx2(globalPos)[dimIdx]) );
        }

        // Production term
        static const auto enableKOmegaProductionLimiter = getParam<bool>("Problem.KOmega.EnableProductionLimiter", false);
        const Scalar wProductionOriginal= 2.0 * nut(globalPos) * shearStressTensorProduct(globalPos);
        Scalar WProductionAlternative = 0.0;
        if (enableKOmegaProductionLimiter)
        {
            // Production limiter (2006 K-limiter via P term)
            WProductionAlternative += 20.0 * turbConstants_.betaK * turbulentKineticEnergy(globalPos) * dissipation(globalPos);
        }
        dissipationSol -= turbConstants_.alpha * std::min(wProductionOriginal, WProductionAlternative) * (dissipation(globalPos) / turbulentKineticEnergy(globalPos));

        // Destruction term
        dissipationSol += turbConstants_.betaOmega * dissipation(globalPos) * dissipation(globalPos);

        // CrossDiffusion term
        static const auto enableCrossDiffusionTerm = getParam<bool>("Problem.KOmega.EnableCrossDiffusionTerm", true);
        if (enableCrossDiffusionTerm)
        {
            Scalar kwProduct = 0.0;
            for (unsigned int j = 0; j < dim; ++j)
                kwProduct +=  dkdx(globalPos)[j] * dwdx(globalPos)[j];

            dissipationSol -= sigma_d(kwProduct) * kwProduct / dissipation(globalPos);
        }

        return dissipationSol;
    }

    //! \brief Returns if the product of \$f  \nabla K \$f and \$f  \nabla omega \$f is more than 0
    bool NabProductBool(const Scalar& NabProduct) const
    { return NabProduct > 0; }

    //! \brief Returns the \$f \sigma_{d} \$f constant
    Scalar sigma_d(const Scalar& NabProduct) const
    { return NabProductBool(NabProduct) == 1 ? 1.0/8.0 : 0; }

    Scalar density_;
    Scalar kinematicViscosity_;
    TurbulentConstants turbConstants_;
};

} // end namespace Dumux

#endif
