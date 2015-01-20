// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Caculates the Jacobian of the local residual for fully-implicit models
 */
#ifndef DUMUX_IMPLICIT_LOCAL_JACOBIAN_HH
#define DUMUX_IMPLICIT_LOCAL_JACOBIAN_HH

#include <dune/common/version.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include "implicitproperties.hh"

namespace Dumux
{
/*!
 * \ingroup ImplicitLocalJacobian
 * \brief Calculates the Jacobian of the local residual for fully-implicit models
 *
 * The default behavior is to use numeric differentiation, i.e.
 * forward or backward differences (2nd order), or central
 * differences (3rd order). The method used is determined by the
 * "NumericDifferenceMethod" property:
 *
 * - if the value of this property is smaller than 0, backward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
 *   \f]
 *
 * - if the value of this property is 0, central
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
 *   \f]
 *
 * - if the value of this property is larger than 0, forward
 *   differences are used, i.e.:
 *   \f[
 \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
 *   \f]
 *
 * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
 * is the value of a sub-control volume's primary variable at the
 * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
 */
template<class TypeTag>
class ImplicitLocalJacobian
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension,

        Green = JacobianAssembler::Green
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::Matrix<MatrixBlock> LocalBlockMatrix;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    // copying a local jacobian is not a good idea
    ImplicitLocalJacobian(const ImplicitLocalJacobian &);

public:
    ImplicitLocalJacobian()
    {
        numericDifferenceMethod_ = GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, NumericDifferenceMethod);
        Valgrind::SetUndefined(problemPtr_);
    }


    /*!
     * \brief Initialize the local Jacobian object.
     *
     * At this point we can assume that everything has been allocated,
     * although some objects may not yet be completely initialized.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;
        localResidual_.init(problem);

        // assume quadrilinears as elements with most vertices
        if (isBox)
        {
            A_.setSize(2<<dim, 2<<dim);
            storageJacobian_.resize(2<<dim);
        }
        else 
        {
            A_.setSize(1, 2<<dim);
            storageJacobian_.resize(1);
        }
    }

    /*!
     * \brief Assemble an element's local Jacobian matrix of the
     *        defect.
     *
     * \param element The DUNE Codim<0> entity which we look at.
     */
    void assemble(const Element &element)
    {
        // set the current grid element and update the element's
        // finite volume geometry
        elemPtr_ = &element;
        fvElemGeom_.update(gridView_(), element);
        reset_();

        bcTypes_.update(problem_(), element_(), fvElemGeom_);

        // set the hints for the volume variables
        model_().setHints(element, prevVolVars_, curVolVars_);

        // update the secondary variables for the element at the last
        // and the current time levels
        prevVolVars_.update(problem_(),
                            element_(),
                            fvElemGeom_,
                            true /* isOldSol? */);

        curVolVars_.update(problem_(),
                           element_(),
                           fvElemGeom_,
                           false /* isOldSol? */);

        // update the hints of the model
        model_().updateCurHints(element, curVolVars_);

        // calculate the local residual
        localResidual().eval(element_(),
                             fvElemGeom_,
                             prevVolVars_,
                             curVolVars_,
                             bcTypes_);
        residual_ = localResidual().residual();
        storageTerm_ = localResidual().storageTerm();

        model_().updatePVWeights(element_(), curVolVars_);

        // calculate the local jacobian matrix
        int numRows, numCols;
        if (isBox)
        {
            numRows = numCols = fvElemGeom_.numScv;
        }
        else 
        {
            numRows = 1;
            numCols = fvElemGeom_.numNeighbors;
        }
        ElementSolutionVector partialDeriv(numRows);
        PrimaryVariables storageDeriv(0.0);
        for (int col = 0; col < numCols; col++) {
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++) {
                asImp_().evalPartialDerivative_(partialDeriv,
                                                storageDeriv,
                                                col,
                                                pvIdx);

                // update the local stiffness matrix with the current partial
                // derivatives
                updateLocalJacobian_(col,
                                     pvIdx,
                                     partialDeriv,
                                     storageDeriv);
            }
        }
    }

    /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    const LocalResidual &localResidual() const
    { return localResidual_; }

    /*!
     * \brief Returns a reference to the object which calculates the
     *        local residual.
     */
    LocalResidual &localResidual()
    { return localResidual_; }

    /*!
     * \brief Returns the Jacobian of the equations at subcontrolvolume i 
     * to the primary variables at subcontrolvolume j.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     * \param j The local subcontrolvolume index which holds
     *          primary variables
     */
    const MatrixBlock &mat(const int i, const int j) const
    { return A_[i][j]; }

    /*!
     * \brief Returns the Jacobian of the storage term at subcontrolvolume i.
     *
     * \param i The local subcontrolvolume index
     */
    const MatrixBlock &storageJacobian(const int i) const
    { return storageJacobian_[i]; }

    /*!
     * \brief Returns the residual of the equations at subcontrolvolume i.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     */
    const PrimaryVariables &residual(const int i) const
    { return residual_[i]; }

    /*!
     * \brief Returns the storage term for subcontrolvolume i.
     *
     * \param i The local subcontrolvolume index on which
     *          the equations are defined
     */
    const PrimaryVariables &storageTerm(const int i) const
    { return storageTerm_[i]; }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param scvIdx     The local index of the element's subcontrolvolume for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon(const int scvIdx,
                          const int pvIdx) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        Scalar priVar = this->curVolVars_[scvIdx].priVar(pvIdx);
        return baseEps*(std::abs(priVar) + 1.0);
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    {
        Valgrind::CheckDefined(problemPtr_);
        return *problemPtr_;
    }

    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the element.
     */
    const Element &element_() const
    {
        Valgrind::CheckDefined(elemPtr_);
        return *elemPtr_;
    }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the jacobian assembler.
     */
    const JacobianAssembler &jacAsm_() const
    { return model_().jacobianAssembler(); }

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Reset the local jacobian matrix to 0
     */
    void reset_()
    {
        for (unsigned int i = 0; i < A_.N(); ++ i) {
            storageJacobian_[i] = 0.0;
            for (unsigned int j = 0; j < A_.M(); ++ j) {
                A_[i][j] = 0.0;
            }
        }
    }

    /*!
     * \brief Compute the partial derivatives to a primary variable at
     *        an degree of freedom.
     *
     * This method can be overwritten by the implementation if a
     * better scheme than numerical differentiation is available.
     *
     * The default implementation of this method uses numeric
     * differentiation, i.e. forward or backward differences (2nd
     * order), or central differences (3rd order). The method used is
     * determined by the "NumericDifferenceMethod" property:
     *
     * - if the value of this property is smaller than 0, backward
     *   differences are used, i.e.:
     *   \f[
         \frac{\partial f(x)}{\partial x} \approx \frac{f(x) - f(x - \epsilon)}{\epsilon}
     *   \f]
     *
     * - if the value of this property is 0, central
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x - \epsilon)}{2 \epsilon}
     *   \f]
     *
     * - if the value of this property is larger than 0, forward
     *   differences are used, i.e.:
     *   \f[
           \frac{\partial f(x)}{\partial x} \approx \frac{f(x + \epsilon) - f(x)}{\epsilon}
     *   \f]
     *
     * Here, \f$ f \f$ is the residual function for all equations, \f$x\f$
     * is the value of a sub-control volume's primary variable at the
     * evaluation point and \f$\epsilon\f$ is a small value larger than 0.
     *
     * \param partialDeriv The vector storing the partial derivatives of all
     *              equations
     * \param storageDeriv the mass matrix contributions
     * \param col The block column index of the degree of freedom 
     *            for which the partial derivative is calculated.
     *            Box: a sub-control volume index.
     *            Cell centered: a neighbor index.
     * \param pvIdx The index of the primary variable 
     *              for which the partial derivative is calculated
     */
    void evalPartialDerivative_(ElementSolutionVector &partialDeriv,
                                PrimaryVariables &storageDeriv,
                                const int col,
                                const int pvIdx)
    {
        int dofIdxGlobal;
        FVElementGeometry neighborFVGeom;
        ElementPointer neighbor(element_());
        if (isBox)
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            dofIdxGlobal = vertexMapper_().subIndex(element_(), col, dim);
#else
            dofIdxGlobal = vertexMapper_().map(element_(), col, dim);
#endif
        }
        else
        {
            neighbor = fvElemGeom_.neighbors[col];
            neighborFVGeom.updateInner(*neighbor);
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            dofIdxGlobal = problemPtr_->elementMapper().index(*neighbor);
#else
            dofIdxGlobal = problemPtr_->elementMapper().map(*neighbor);
#endif
        }

        PrimaryVariables priVars(model_().curSol()[dofIdxGlobal]);
        VolumeVariables origVolVars(curVolVars_[col]);

        curVolVars_[col].setEvalPoint(&origVolVars);
        Scalar eps = asImp_().numericEpsilon(col, pvIdx);
        Scalar delta = 0;

        if (numericDifferenceMethod_ >= 0) {
            // we are not using backward differences, i.e. we need to
            // calculate f(x + \epsilon)

            // deflect primary variables
            priVars[pvIdx] += eps;
            delta += eps;

            // calculate the residual
            if (isBox)
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        element_(),
                                        fvElemGeom_,
                                        col,
                                        false);
            else
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        *neighbor,
                                        neighborFVGeom,
                                        /*scvIdx=*/0,
                                        false);
                   
            localResidual().eval(element_(),
                                 fvElemGeom_,
                                 prevVolVars_,
                                 curVolVars_,
                                 bcTypes_);

            // store the residual and the storage term
            partialDeriv = localResidual().residual();
            if (isBox || col == 0)
                storageDeriv = localResidual().storageTerm()[col];
        }
        else {
            // we are using backward differences, i.e. we don't need
            // to calculate f(x + \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv = residual_;
            if (isBox || col == 0)
                storageDeriv = storageTerm_[col];
        }


        if (numericDifferenceMethod_ <= 0) {
            // we are not using forward differences, i.e. we
            // need to calculate f(x - \epsilon)

            // deflect the primary variables
            priVars[pvIdx] -= delta + eps;
            delta += eps;

            // calculate residual again
            if (isBox)
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        element_(),
                                        fvElemGeom_,
                                        col,
                                        false);
            else
                curVolVars_[col].update(priVars,
                                        problem_(),
                                        *neighbor,
                                        neighborFVGeom,
                                        /*scvIdx=*/0,
                                        false);

            localResidual().eval(element_(),
                                 fvElemGeom_,
                                 prevVolVars_,
                                 curVolVars_,
                                 bcTypes_);
            partialDeriv -= localResidual().residual();
            if (isBox || col == 0)
                storageDeriv -= localResidual().storageTerm()[col];
        }
        else {
            // we are using forward differences, i.e. we don't need to
            // calculate f(x - \epsilon) and we can recycle the
            // (already calculated) residual f(x)
            partialDeriv -= residual_;
            if (isBox || col == 0)
                storageDeriv -= storageTerm_[col];
        }

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        partialDeriv /= delta;
        storageDeriv /= delta;

        // restore the original state of the element's volume variables
        curVolVars_[col] = origVolVars;

#if HAVE_VALGRIND
        for (unsigned i = 0; i < partialDeriv.size(); ++i)
            Valgrind::CheckDefined(partialDeriv[i]);
#endif
    }

    /*!
     * \brief Updates the current local Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col' .
     */
    void updateLocalJacobian_(const int col,
                              const int pvIdx,
                              const ElementSolutionVector &partialDeriv,
                              const PrimaryVariables &storageDeriv)
    {
        // store the derivative of the storage term
        if (isBox || col == 0)
        {
            for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                storageJacobian_[col][eqIdx][pvIdx] = storageDeriv[eqIdx];
            }
        }

        for (int i = 0; i < fvElemGeom_.numScv; i++)
        {
            // Green vertices are not to be changed!
            if (!isBox || jacAsm_().vertexColor(element_(), i) != Green) {
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++) {
                    // A[i][col][eqIdx][pvIdx] is the rate of change of
                    // the residual of equation 'eqIdx' at dof 'i'
                    // depending on the primary variable 'pvIdx' at dof
                    // 'col'.
                    this->A_[i][col][eqIdx][pvIdx] = partialDeriv[i][eqIdx];
                    Valgrind::CheckDefined(this->A_[i][col][eqIdx][pvIdx]);
                }
            }
        }
    }

    const Element *elemPtr_;
    FVElementGeometry fvElemGeom_;

    ElementBoundaryTypes bcTypes_;

    // The problem we would like to solve
    Problem *problemPtr_;

    // secondary variables at the previous and at the current time
    // levels
    ElementVolumeVariables prevVolVars_;
    ElementVolumeVariables curVolVars_;

    LocalResidual localResidual_;

    LocalBlockMatrix A_;
    std::vector<MatrixBlock> storageJacobian_;

    ElementSolutionVector residual_;
    ElementSolutionVector storageTerm_;

    int numericDifferenceMethod_;
};
}

#endif
