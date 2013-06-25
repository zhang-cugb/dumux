//$Id: el2pnewtoncontroller.hh 10663 2013-05-14 16:22:40Z melanie $
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief An el2p specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_EL2P_NEWTON_CONTROLLER_HH
#define DUMUX_EL2P_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {

template <class TypeTag>
class ElTwoPNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;

    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

public:
    /*!
     * \brief Destructor
     */
    ElTwoPNewtonController(const Problem &problem)
    : ParentType(problem),linearSolver_(problem)
    {
        this->setTargetSteps(9);
        this->setMaxSteps(18);
    };

    void newtonUpdateRelError(const SolutionVector &uOld,
                              const SolutionVector &deltaU)
    {
        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        this->error_ = 0;

        for (int i = 0; i < int(uOld.size()); ++i) {
            Scalar vertErr = std::abs(deltaU[i]/(1.0 + std::abs((uOld[i]) + uOld[i] - deltaU[i])/2));
            this->error_ = std::max(this->error_, vertErr);
        }

        this->error_ = this->gridView_().comm().max(this->error_);
    }

    void newtonUpdate(SolutionVector &uCurrentIter,
            const SolutionVector &uLastIter,
            const SolutionVector &deltaU)
    {
//        this->writeConvergence_(uLastIter, deltaU);

        newtonUpdateRelError(uLastIter, deltaU);

        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;

//        printvector(std::cout, deltaU, "new solution", "row", 12, 1, 3);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    void newtonSolveLinear(const JacobianMatrix &A,
                           SolutionVector &x,
                           const SolutionVector &b)
    {
        try {
            if (this->numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                if (this->gridView_().comm().size() > 1)
                    norm2 = this->gridView_().comm().sum(norm2);

                initialAbsoluteError_ = std::sqrt(norm2);
                lastAbsoluteError_ = initialAbsoluteError_;
            }

            int converged = linearSolver_.solve(A, x, b);
//            printvector(std::cout, x.base(), "x", "row", 5, 1, 5);
//            printvector(std::cout, b.base(), "rhs", "row", 5, 1, 5);
//            Dune::writeMatrixToMatlab(A.base(), "matrix.txt");

            // make sure all processes converged
            int convergedRemote = converged;
            if (this->gridView_().comm().size() > 1)
                convergedRemote = this->gridView_().comm().min(converged);

            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge");
            }
            else if (!convergedRemote) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge on a remote process");
            }
        }
        catch (Dune::MatrixBlockError e) {
            // make sure all processes converged
            int converged = 0;
            if (this->gridView_().comm().size() > 1)
                converged = this->gridView_().comm().min(converged);

            Dumux::NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (this->gridView_().comm().size() > 1)
                converged = this->gridView_().comm().min(converged);

            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    // absolute errors and tolerance
    Scalar absoluteError_;
    Scalar lastAbsoluteError_;
    Scalar initialAbsoluteError_;
    Scalar absoluteTolerance_;

    // the linear solver
    LinearSolver linearSolver_;

};
}

#endif
