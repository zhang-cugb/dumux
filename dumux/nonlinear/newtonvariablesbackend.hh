// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Nonlinear
 * \brief Backends for the Newton for different variable types
 */
#ifndef DUMUX_NEWTON_VARIABLES_BACKEND_HH
#define DUMUX_NEWTON_VARIABLES_BACKEND_HH

#include <type_traits>
#include <dune/common/typetraits.hh>
#include <dune/common/std/type_traits.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief Class providing operations with variables needed in the Newton solver
 */
template<class DofVector, class Enable = void>
class NewtonDofBackend;

/*!
 * \ingroup Nonlinear
 * \brief Class providing operations with variables needed in the Newton solver
 */
// template<class Variables, class Enable = void>
// class NewtonVariablesBackend;

/*!
 * \file
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 */
template<class Scalar>
class NewtonDofBackend<Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, Scalar>>
{
public:
    using DofVector = Scalar; //!< the type of the dofs parametrizing the variables object

    static std::size_t size(const DofVector& d)
    { return 1; }

    static DofVector makeDofVector(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVector(std::size_t size)
    { return 0.0; }

    static Scalar makeDofVectorForSolver(const DofVector& d)
    { return d; }

    static Scalar makeZeroDofVectorForSolver(std::size_t size)
    { return 0.0; }

    static Scalar reconstructDofVectorFromSolver(const Scalar& d)
    { return d; }

    static Scalar maxRelativeShift(const DofVector& previous, const DofVector& delta)
    {
        const auto current = previous - delta;
        using std::abs; using std::max;
        auto error = abs(previous - current);
        error /= max<Scalar>(1.0, abs(previous + current)/2.0);
        return error;
    }
};

/*!
 * \file
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 */
template<class BT>
class NewtonDofBackend<Dune::BlockVector<BT>>
{
    // TODO: CHECK PRIVARSWITCH VARIABLES

    using Scalar = typename BT::value_type;

public:
    using DofVector = Dune::BlockVector<BT>; //!< the type of the dofs parametrizing the variables object

    static std::size_t size(const DofVector& d)
    { return d.size(); }

    static DofVector makeDofVector(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVector(std::size_t size)
    { DofVector d; d.resize(size); return d; }

    static DofVector makeDofVectorForSolver(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVectorForSolver(std::size_t size)
    { DofVector d; d.resize(size); return d; }

    static DofVector reconstructDofVectorFromSolver(const DofVector& d)
    { return d; }

    static Scalar maxRelativeShift(const DofVector& previous, const DofVector& delta)
    {
        Scalar shift = 0.0;
        for (int i = 0; i < int(previous.size()); ++i)
        {
            auto uNew = previous[i];
            uNew -= delta[i];

            // TODO: impl backend for FieldVector
            for (unsigned int j = 0; j < uNew.size(); ++j)
            {
                using std::max;
                using std::abs;
                auto error = abs(previous[i][j] - uNew[j]);
                error /= max<Scalar>(1.0, abs(previous[i][j] + uNew[j])/2.0);
                shift = max(shift, error);

                // using ScalarBackend = NewtonDofBackend<Scalar, Scalar>;
                // shift = max(shift, ScalarBackend::maxRelativeShift(previous[i][j], uNew[j]));
            }
        }

        return shift;
    }
};

namespace Impl {

template<class Vars>
using SolutionVectorType = typename Vars::SolutionVector;

template<class Vars, bool HasSolVec>
class NewtonVariablesBackend;

/*!
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 */
template<class Vars>
class NewtonVariablesBackend<Vars, false>
: public NewtonDofBackend<Vars>
{
    using ParentType = NewtonDofBackend<Vars>;

public:
    using Variables = Vars;
    using typename ParentType::DofVector;

    static void update(Variables& v, const DofVector& dofs)
    { v = dofs; }

    //! operations on variables
    static const DofVector& getDofVector(Variables& v)
    { return v; }
};

/*!
 * \file
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 */
template<class Vars>
class NewtonVariablesBackend<Vars, true>
: public NewtonDofBackend<typename Vars::SolutionVector>
{
public:
    using DofVector = typename Vars::SolutionVector;
    using Variables = Vars; //!< the type of the variables object

    static void update(Variables& v, const DofVector& dofs)
    { v.updateDofs(dofs); }

    //! operations on variables
    static const DofVector& getDofVector(Variables& v)
    { return v.dofs(); }
};
} // end namespace Impl

/*!
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 */
template<class Vars>
using NewtonVariablesBackend = Impl::NewtonVariablesBackend<Vars, Dune::Std::is_detected_v<Impl::SolutionVectorType, Vars>>;

} // end namespace Dumux

#endif
