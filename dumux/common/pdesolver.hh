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
 * \ingroup Common
 * \brief Defines a high-level interface for a PDESolver
 */
#ifndef DUMUX_COMMON_PDESOLVER_HH
#define DUMUX_COMMON_PDESOLVER_HH

#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>

// forward declare
namespace Dune {
template <class FirstRow, class ... Args>
class MultiTypeBlockMatrix;
}

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A high-level interface for a PDESolver
 *
 * A PDESolver is constructed with an assembler and a linear solver
 * and has a method solve that linearizes (if not already linear), assembles, solves and updates
 * given an initial solution producing a new solution.
 *
 * \tparam AssemblerType A PDE linearized system assembler
 * \tparam LinearSolverType A linear system solver
 */
template<class AssemblerType, class LinearSolverType>
class PDESolver
{
    using SolutionVector = typename AssemblerType::ResidualType;
    using Variables = typename AssemblerType::Variables;
    using Scalar = typename AssemblerType::Scalar;

public:
    //! export the underlying types for assembler and linear solver
    using Assembler = AssemblerType;
    using LinearSolver = LinearSolverType;

    //! The constructor
    PDESolver(std::shared_ptr<Assembler> assembler,
              std::shared_ptr<LinearSolver> linearSolver)
    : assembler_(assembler)
    , linearSolver_(linearSolver)
    {}

    virtual ~PDESolver() = default;

    /*!
     * \brief Solve the given PDE system (usually assemble + solve linear system + update)
     * \param sol a solution vector possbilty containing an initial solution
     * \param vars an instance of the variables object containing all quantities
     *             required to evaluate the PDE (may be more than just the primary variables)
     * \note The variables depend on the primary variables. Thus, here we assume that the
     *       provided solution is the one which has been used to construct the variables object.
     */
    virtual void solve(SolutionVector& sol, Variables& vars) = 0;

    /*!
     * \brief Access the assembler
     */
    const Assembler& assembler() const
    { return *assembler_; }

    /*!
     * \brief Access the assembler
     */
    Assembler& assembler()
    { return *assembler_; }

    /*!
     * \brief Access the linear solver
     */
    const LinearSolver& linearSolver() const
    { return *linearSolver_; }

    /*!
     * \brief Access the linear solver
     */
    LinearSolver& linearSolver()
    { return *linearSolver_; }

protected:

    /*!
     * \brief Helper function to assure the MultiTypeBlockMatrix's sub-blocks have the correct sizes.
     */
    template <class FirstRow, class ... Args>
    bool checkSizesOfSubMatrices(const Dune::MultiTypeBlockMatrix<FirstRow, Args...>& matrix) const
    {
        bool matrixHasCorrectSize = true;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Dune::MultiTypeBlockMatrix<FirstRow, Args...>::N()>(), [&](const auto i)
        {
            const auto& row = matrix[i];
            const auto numRowsLeftMostBlock = row[Dune::index_constant<0>{}].N();
            forEach(row, [&](const auto& subBlock)
            {
                if (subBlock.N() != numRowsLeftMostBlock)
                    matrixHasCorrectSize = false;
            });
        });
        return matrixHasCorrectSize;
    }

private:
    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<LinearSolver> linearSolver_;
};

} // namespace Dumux

#endif
