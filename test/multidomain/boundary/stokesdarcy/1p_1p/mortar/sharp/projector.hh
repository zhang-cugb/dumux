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
#ifndef DUMUX_MORTAR_SHARP_PROJECTOR_HH
#define DUMUX_MORTAR_SHARP_PROJECTOR_HH

#include <type_traits>

#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/io.hh>

#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux {

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class SolutionVector>
class SharpMortarProjector
{
    using Scalar = typename SolutionVector::field_type;
    using BlockType = typename SolutionVector::block_type;
    using SubMatrix = Dune::BCRSMatrix< Dune::FieldMatrix<Scalar, 1, 1> >;
    static_assert(BlockType::size() == 1, "Currently only block size of 1 is supported");

    // Helper aliases
    template<class... ColTypes> using MatrixRow = Dune::MultiTypeBlockVector<ColTypes...>;
    template<class... RowTypes> using Matrix = Dune::MultiTypeBlockMatrix<RowTypes...>;

    // Aliases for the matrix and vector types
    using Vector2 = Dune::MultiTypeBlockVector<SolutionVector, SolutionVector>;
    using Vector3 = Dune::MultiTypeBlockVector<SolutionVector, SolutionVector, SolutionVector>;

    using MatrixRow2 = MatrixRow<SubMatrix, SubMatrix>;
    using MatrixRow3 = MatrixRow<SubMatrix, SubMatrix, SubMatrix>;

    using Matrix23 = Matrix<MatrixRow3, MatrixRow3>;
    using Matrix33 = Matrix<MatrixRow3, MatrixRow3, MatrixRow3>;
    using Matrix32 = Matrix<MatrixRow2, MatrixRow2, MatrixRow2>;

public:
    /*!
     * \brief TODO doc me.
     */
    template<class Projector>
    SharpMortarProjector(const Projector& toSubDomain1,
                         const Projector& toSubDomain2)
    {
        const auto& M1 = toSubDomain1.massMatrix();
        const auto& M2 = toSubDomain2.massMatrix();

        const auto& B1 = toSubDomain1.projectionMatrix();
        const auto& B2 = toSubDomain2.projectionMatrix();

        // create identity matrices
        SubMatrix I1, I2;

        Dune::MatrixIndexSet pattern1, pattern2;
        pattern1.resize(B1.N(), B1.N());
        pattern2.resize(B2.N(), B2.N());
        for (std::size_t i = 0; i < B1.N(); ++i)
            pattern1.add(i, i);
        for (std::size_t i = 0; i < B2.N(); ++i)
            pattern2.add(i, i);

        pattern1.exportIdx(I1); I1 = 1.0;
        pattern2.exportIdx(I2); I2 = 1.0;

        // set up the transposed matrices using identities
        SubMatrix B1T, B2T;
        Dune::transposeMatMultMat(B1T, B1, I1);
        Dune::transposeMatMultMat(B2T, B2, I2);

        using namespace Dune::Indices;

        // fill the sub matrices
        I_[_0][_0] = I1;          I_[_0][_1] = SubMatrix{}; I_[_0][_2] = SubMatrix{};
        I_[_1][_0] = SubMatrix{}; I_[_1][_1] = I2;          I_[_1][_2] = SubMatrix{};

        A_[_0][_0] = M1;          A_[_0][_1] = SubMatrix{}; A_[_0][_2] = B1;
        A_[_1][_0] = SubMatrix{}; A_[_1][_1] = M2;          A_[_1][_2] = B2;
        A_[_2][_0] = B1T;         A_[_2][_1] = B2T;         A_[_2][_2] = SubMatrix{};

        C_[_0][_0] = B1;          C_[_0][_1] = SubMatrix{};
        C_[_1][_0] = SubMatrix{}; C_[_1][_1] = B2;
        C_[_2][_0] = SubMatrix{}; C_[_2][_1] = SubMatrix{};

        // set up the zero sub-matrix block sizes
        using IndexSet = Dune::MatrixIndexSet;

        IndexSet{I1.N(), I2.M()}.exportIdx(I_[_0][_1]); // I[0][1]
        IndexSet{I1.N(), B1.M()}.exportIdx(I_[_0][_2]); // I[0][2]
        IndexSet{I2.N(), I1.M()}.exportIdx(I_[_1][_0]); // I[1][0]
        IndexSet{I2.N(), B2.M()}.exportIdx(I_[_1][_2]); // I[1][2]

        IndexSet{M1.N(), M2.M()}.exportIdx(A_[_0][_1]); // A[0][1]
        IndexSet{M2.N(), M1.M()}.exportIdx(A_[_1][_0]); // A[1][0]
        IndexSet{B1.M(), B1.M()}.exportIdx(A_[_2][_2]); // A[2][2]

        IndexSet{B1.N(), B2.M()}.exportIdx(C_[_0][_1]); // C[0][1]
        IndexSet{B2.N(), B1.M()}.exportIdx(C_[_1][_0]); // C[1][0]
        IndexSet{B1.M(), B1.M()}.exportIdx(C_[_2][_0]); // C[2][0]
        IndexSet{B2.M(), B2.M()}.exportIdx(C_[_2][_1]); // C[2][1]
    }

    /*!
     * \brief TODO doc me.
     */
    Vector2 projectMortarToSubDomain(const SolutionVector& x1,
                                     const SolutionVector& x2) const
    {
        using namespace Dune::Indices;

        Vector2 x(x1, x2);

        // set up right hand side rhs = Cx
        Vector3 rhs;
        rhs[_0].resize(C_[_0][_0].N());
        rhs[_1].resize(C_[_1][_0].N());
        rhs[_2].resize(C_[_2][_0].N());
        C_.mv(x, rhs);

        // solve A*xTmp = rhs
        // solver doesn't support MultyTypeMatrices,
        // so we copy the values beforehand
        auto ATmp = MatrixConverter<Matrix33>::multiTypeToBCRSMatrix(A_);
        auto rhsTmp = VectorConverter<Vector3>::multiTypeToBlockVector(rhs);
        auto tmp = rhsTmp;

        UMFPackBackend solver;
        solver.solve(ATmp, tmp, rhsTmp);

        // copy values back into MultiTypeVector (reuse rhs)
        VectorConverter<Vector3>::retrieveValues(rhs, tmp);

        // restrict
        Vector2 xp;
        xp[_0].resize(I_[_0][_0].N());
        xp[_1].resize(I_[_1][_0].N());
        I_.mv(rhs, xp);

        return xp;
    }

    /*!
     * \brief TODO doc me.
     */
    Vector2 projectSubDomainToMortar(const SolutionVector& x1,
                                     const SolutionVector& x2) const
    {
        using namespace Dune::Indices;

        // set up right hand side
        Vector2 x(x1, x2);

        Vector3 rhs;
        rhs[_0].resize(I_[_0][_0].M());
        rhs[_1].resize(I_[_0][_1].M());
        rhs[_2].resize(I_[_0][_2].M());

        // compute rhs = I_^T*x
        I_[_0][_0].mtv(x[_0], rhs[_0]);
        I_[_1][_0].umtv(x[_1], rhs[_0]);
        I_[_0][_1].mtv(x[_0], rhs[_1]);
        I_[_1][_1].umtv(x[_1], rhs[_1]);
        I_[_0][_2].mtv(x[_0], rhs[_2]);
        I_[_1][_2].umtv(x[_1], rhs[_2]);

        // // solve A*xTmp = rhs
        // // solver doesn't support MultyTypeMatrices,
        // // so we copy the values beforehand
        // auto ATmp = MatrixConverter<Matrix33>::multiTypeToBCRSMatrix(A_);
        // auto rhsTmp = VectorConverter<Vector3>::multiTypeToBlockVector(rhs);
        // auto tmp = rhsTmp;
        //
        // UMFPackBackend solver;
        // solver.solve(ATmp, tmp, rhsTmp);
        //
        // // copy values back into MultiTypeVector (reuse rhs)
        // VectorConverter<Vector3>::retrieveValues(rhs, tmp);

        // apply C
        Vector2 xp;
        xp[_0].resize(C_[_0][_0].M());
        xp[_1].resize(C_[_0][_1].M());

        C_[_0][_0].mtv(rhs[_0], xp[_0]);
        C_[_1][_0].umtv(rhs[_1], xp[_0]);
        C_[_2][_0].umtv(rhs[_2], xp[_0]);

        C_[_0][_1].mtv(rhs[_0], xp[_1]);
        C_[_1][_1].umtv(rhs[_1], xp[_1]);
        C_[_2][_1].umtv(rhs[_2], xp[_1]);

        return xp;
    }

private:
    Matrix23 I_;
    Matrix33 A_;
    Matrix32 C_;
};

} // end namespace Dumux

#endif
