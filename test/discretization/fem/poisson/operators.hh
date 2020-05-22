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
 * \todo TODO: WHICH GROUP?
 * \ingroup TODO: WHICH GROUP?
 * \copydoc FEPoissonOperators
 */
#ifndef DUMUX_FEM_POISSON_OPERATORS_HH
#define DUMUX_FEM_POISSON_OPERATORS_HH

#include <dumux/common/math.hh>
#include <dumux/assembly/fem/operatorsbase.hh>

namespace Dumux {

/*!
 * \file
 * \todo TODO: WHICH GROUP?
 * \ingroup TODO: WHICH GROUP?
 * \brief The operators class for a poisson problem
 *        using the finite element method.
 * \tparam GridVarsLocalView The type of local view on the grid variables
 */
template<class GridVarsLocalView>
class FEPoissonOperators
: public FEOperatorsBase< GridVarsLocalView >
{
    using ParentType = FEOperatorsBase<GridVarsLocalView>;
    using IpVariables = typename GridVarsLocalView::GridVariables::IntegrationPointVariables;

public:
    //! export flux term type
    using typename ParentType::FluxTerm;

    //! pull up base class constructors
    using ParentType::ParentType;

    /*!
     * \brief Calculate the flux term of the equation
     * \param elemSol The element solution vector
     * \param ipData The shape function values/gradients evaluated at the integration point
     * \param ipVars The primary/secondary variables evaluated at the integration point
     */
    template<class IpData>
    FluxTerm flux(const IpData& ipData, const IpVariables& ipVars) const
    {
        // evaluate gradient in solution
        typename IpData::GlobalPosition gradX(0.0);
        for (unsigned int i = 0; i < ipData.size(); ++i)
        {
            auto tmp = ipData.gradN(i);
            tmp *= this->gridVariablesLocalView().elemSol()[i];
            gradX += tmp;
        }

        // The flux is tensor*gradX
        FluxTerm result(0.0);
        assert(result.size() == 1);
        result[0] = mv(this->problem_().poissonTensor(), gradX);

        return result;
    }
};

} // end namespace Dumux

#endif
