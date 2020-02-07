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
 * \ingroup AnalyticalRans
 * \brief
 */
#ifndef DUMUX_FREEFLOW_RANS_KOMEGA_ANALYTICAL_HH
#define DUMUX_FREEFLOW_RANS_KOMEGA_ANALYTICAL_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace KOmegaAnalytical {

template<class Scalar, class GlobalPostion, int dim>
class KOmegaAnalytical
{
public:
    using VelocityVector = std::array<Scalar, dim>;
    using VelocityGradientMatrix = std::array<VelocityVector, dim>;
    using ScalarGradientVector = std::array<Scalar, dim>;

    static VelocityVector velocity(const GlobalPostion& globalPos) const
    {
        VelocityVector v(0.0);
        v[0] = 2.0 * globalPos[0] * globalPos[0] * globalPos[0];
        if (dim > 1)
            v[1] = 2.0 - 2.0 * globalPos[1] * globalPos[1] * globalPos[1];
        return v;
    }

    static Scalar pressure(const GlobalPostion& globalPos) const
    {
        if (dim > 1)
            return 2.0 - 2.0 * globalPos[0] * globalPos[1];
        else
            return 2.0 - 2.0 * globalPos[0];
    }

    static Scalar turbulentKineticEnergy(const GlobalPostion& globalPos) const
    {
        if (dim > 1)
            return 1.0 + globalPos[0] * globalPos[0] * globalPos[1] * globalPos[1];
        else
            return 1.0 + globalPos[0] * globalPos[0];
    }

    static Scalar dissipation(const GlobalPostion& globalPos) const
    {
        if (dim > 1)
            return 2.0 - globalPos[0] * globalPos[0] * globalPos[1] * globalPos[1];
        else
            return 2.0 - globalPos[0] * globalPos[0];
    }

    //! \brief The gradients
    static VelocityGradientMatrix dvdx(const GlobalPostion& globalPos) const
    {
        VelocityGradientMatrix dvdx;
        dvdx[0][0] = 6.0 * globalPos[0] * globalPos[0];
        dvdx[0][1] = 0.0;
        dvdx[1][0] = 0.0;
        dvdx[1][1] = -6.0 * globalPos[1] * globalPos[1];

        return dvdx;
    }

    static VelocityGradientMatrix dv2dx(const GlobalPostion& globalPos) const
    {
        VelocityGradientMatrix dv2dx;
        for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                dv2dx[velIdx][dimIdx] = dvdx(globalPos)[velIdx][dimIdx] * velocity(globalPos)[dimIdx]
                                      + dvdx(globalPos)[dimIdx][dimIdx] * velocity(globalPos)[velIdx];
        return dv2dx;
    }

    static VelocityGradientMatrix dvdx2(const GlobalPostion& globalPos) const
    {
        VelocityGradientMatrix dvdx2;
        dvdx2[0][0] = 12.0 * globalPos[0];
        if (dim > 1)
        {
            dvdx2[0][1] = 0.0;
            dvdx2[1][0] = 0.0;
            dvdx2[1][1] = -12.0 * globalPos[1];
        }
        return dvdx2;
    }

    static ScalarGradientVector dpdx(const GlobalPostion& globalPos) const
    {
        ScalarGradientVector dpdx(0.0);
        dpdx[0] = -2.0;
        if (dim > 1)
        {
            dpdx[0] = -2.0 * globalPos[1];
            dpdx[1] = -2.0 * globalPos[0];
        }
        return dpdx;
    }

    static ScalarGradientVector dkdx(const GlobalPostion& globalPos) const
    {
        ScalarGradientVector dkdx(0.0);
        dkdx[0] = 2.0 * globalPos[0];
        if (dim > 1)
        {
            dkdx[0] = 2.0 * globalPos[0] * globalPos[1] * globalPos[1];
            dkdx[1] = 2.0 * globalPos[1] * globalPos[0] * globalPos[0];
        }
        return dkdx;
    }

    static ScalarGradientVector dkdx2(const GlobalPostion& globalPos) const
    {
        ScalarGradientVector dkdx2(0.0);
        dkdx2[0] = 2.0;
        if (dim > 1)
        {
            dkdx2[0] = 2.0 * globalPos[1] * globalPos[1];
            dkdx2[1] = 2.0 * globalPos[0] * globalPos[0];
        }
        return dkdx2;
    }

    static ScalarGradientVector dwdx(const GlobalPostion& globalPos) const
    {
        ScalarGradientVector dwdx(0.0);
        dwdx[0] = -2.0 * globalPos[0];
        if (dim > 1)
        {
            dwdx[0] = -2.0 * globalPos[0] * globalPos[1] * globalPos[1];
            dwdx[1] = -2.0 * globalPos[0] * globalPos[0] * globalPos[1];
        }
        return dwdx;
    }

    static ScalarGradientVector dwdx2(const GlobalPostion& globalPos) const
    {
        ScalarGradientVector dwdx2(0.0);
        dwdx2[0] = -2.0;
        if (dim > 1)
        {
            dwdx2[0] = -2.0 * globalPos[1] * globalPos[1];
            dwdx2[1] = -2.0 * globalPos[0] * globalPos[0];
        }
        return dwdx2;
    }

    // eddy viscosity and its gradient
    static Scalar nut(const GlobalPostion& globalPos) const
    { return turbulentKineticEnergy(globalPos) / dissipation(globalPos); }

    static ScalarGradientVector dnutdx(const GlobalPostion& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        ScalarGradientVector dnutdx(0.0);
        dnutdx[0] = 6*x
                  / ((2-(x*x))
                  * (2-(x*x)));
        if (dim > 1)
        {
            dnutdx[0] = 6*x*y*y
                      / ((2-(x*x*y*y))
                      * (2-(x*x*y*y)));
            dnutdx[1] = 6*x*x*y
                      / ((2-(x*x*y*y))
                      * (2-(x*x*y*y)));
        }
        return dnutdx;
    }

    // shearStressTensorProduct
    static Scalar shearStressTensorProduct(const GlobalPostion& globalPos) const
    {
        Scalar shearStressTensorProduct = 0.0;
        for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                shearStressTensorProduct += 0.5 * (dvdx(globalPos)[velIdx][dimIdx] + dvdx(globalPos)[dimIdx][velIdx])
                                          * 0.5 * (dvdx(globalPos)[velIdx][dimIdx] + dvdx(globalPos)[dimIdx][velIdx]);
        return shearStressTensorProduct;
    }

    //! \brief Source term values
    static Scalar sourceMassBalanceAtPos(const GlobalPostion& globalPos) const
    {
        // term div(rho*v)
        Scalar v(0.0);
        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            v+= dvdx(global)[dimIdx][dimIdx];

        return density_ * v;
    }

    VelocityVector sourceMomentumBalanceAtPos(const GlobalPostion& globalPos) const
    {
        DimVector y(0.0);
        DimVector gravity(0.0);
        gravity[dim-1] = -9.81;

        for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
        {
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                // inertia term
                y[velIdx] += density_ * dv2dx(global)[velIdx][dimIdx];

                // viscous term (molecular)
                y[velIdx] -= density_ * kinematicViscosity_* dvdx2(global)[velIdx][dimIdx];
                y[velIdx] -= density_ * kinematicViscosity_* dvdx2(global)[dimIdx][velIdx];

                // viscous term (turbulent) (using product rule)
                y[velIdx] -= density_ * (dnutdx(global)[velIdx] * dvdx(global)[velIdx][dimIdx]
                                          + nut(global) * dvdx2(global)[velIdx][dimIdx]);
                y[velIdx] -= density_ * (dnutdx(global)[dimIdx] * dvdx(global)[velIdx][dimIdx]
                                            + nut(global) * dvdx2(global)[dimIdx][velIdx]);
            }

            // pressure term
            y[velIdx] += dpdx(global)[velIdx];

            // gravity term
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
                y[velIdx] -= density_ * gravity[velIdx];
        }
        return y;
    }

};

} // end namespace KOmegaAnalytical
} // end namespace Dumux

#endif
