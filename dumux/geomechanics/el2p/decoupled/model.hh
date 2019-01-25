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
*
* \brief Adaption of the fully implicit scheme to the two-phase linear elasticity model.
*/

#ifndef DUMUX_DECOUPLED_ELASTIC_MODEL_HH
#define DUMUX_DECOUPLED_ELASTIC_MODEL_HH

#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dumux/common/eigenvalues.hh>
#include "dumux/geomechanics/el2p/properties.hh"

namespace Dumux {

namespace Properties {
NEW_PROP_TAG(InitialDisplacement); //!< The initial displacement function
NEW_PROP_TAG(InitialPressSat); //!< The initial pressure and saturation function
}

/*!
 * \ingroup DecoupledElasticModel
 * \brief Adaption of the fully implicit scheme to the two-phase linear elasticity model.
 *
 * This model implements a two-phase flow of compressible immiscible fluids \f$\alpha \in \{ w, n \}\f$.
 * The deformation of the solid matrix is described with a quasi-stationary momentum balance equation.
 * The influence of the pore fluid is accounted for through the effective stress concept (Biot 1941).
 * The total stress acting on a rock is partially supported by the rock matrix and partially supported
 * by the pore fluid. The effective stress represents the share of the total stress which is supported
 * by the solid rock matrix and can be determined as a function of the strain according to Hooke's law.
 *
 * As an equation for the conservation of momentum within the fluid phases the standard multiphase Darcy's approach is used:
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
\f[
 \frac{\partial \phi_{eff} \varrho_\alpha S_\alpha}{\partial t}
 - \text{div} \left\{ \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}
 \mathbf{K}_\text{eff} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 - \phi_{eff} \varrho_\alpha S_\alpha \frac{\partial \mathbf{u}}{\partial t}
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 *
 * A quasi-stationary momentum balance equation is solved for the changes with respect to the initial conditions (Darcis 2012), note
 * that this implementation assumes the soil mechanics sign convention (i.e. compressive stresses are negative):
 \f[
 \text{div}\left( \boldsymbol{\Delta \sigma'}- \Delta p_{eff} \boldsymbol{I} \right) + \Delta \varrho_b {\textbf g} = 0 \;,
 \f]
 * with the effective stress:
 \f[
  \boldsymbol{\sigma'} = 2\,G\,\boldsymbol{\epsilon} + \lambda \,\text{tr} (\boldsymbol{\epsilon}) \, \mathbf{I}.
 \f]
 *
 * and the strain tensor \f$\boldsymbol{\epsilon}\f$ as a function of the solid displacement gradient \f$\textbf{grad} \mathbf{u}\f$:
 \f[
  \boldsymbol{\epsilon} = \frac{1}{2} \, (\textbf{grad} \mathbf{u} + \textbf{grad}^T \mathbf{u}).
 \f]
 *
 * Here, the rock mechanics sign convention is switch off which means compressive stresses are < 0 and tensile stresses are > 0.
 * The rock mechanics sign convention can be switched on for the vtk output via the property system.
 *
 * The effective porosity and the effective permeability are calculated as a function of the solid displacement:
 \f[
      \phi_{eff} = \frac{\phi_{init} + \text{div} \mathbf{u}}{1 + \text{div} \mathbf{u}}
 \f]
 \f[
      K_{eff} = K_{init} \text{exp}\left( 22.2(\phi_{eff}/\phi_{init} -1 )\right)
 \f]
 * The mass balance equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial and the implicit Euler method as time discretization.
 * The momentum balance equations are discretized using a standard Galerkin Finite Element method as
 * spatial discretization scheme.
 *
 *
 * The primary variables are the wetting phase pressure \f$p_w\f$, the nonwetting phase saturation \f$S_n\f$ and the solid
 * displacement vector \f$\mathbf{u}\f$ (changes in solid displacement with respect to initial conditions).
 */
template<class TypeTag>
class DecoupledElasticModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename Element::Geometry::JacobianInverseTransposed JacobianInverseTransposed;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
    typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

public:

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param entity The Entity
     *
     * Due to the mixed discretization schemes which are combined via pdelab for this model
     * the solution vector has a different form than in the pure box models
     * it sorts the primary variables in the following way:
     * p_vertex0 S_vertex0 p_vertex1 S_vertex1 p_vertex2 ....p_vertexN S_vertexN
     * ux_vertex0 uy_vertex0 uz_vertex0 ux_vertex1 uy_vertex1 uz_vertex1 ...
     *
     * Therefore, the serializeEntity function has to be modified.
     */
    template <class Entity>
    void serializeEntity(std::ostream &outStream,
                         const Entity &entity)
    {
        // vertex index
        int dofIdxGlobal = this->dofMapper().index(entity);

        // write phase state
        if (!outStream.good()) {
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << dofIdxGlobal);
        }
        int numScv = this->gridView().size(dim);
        // get p and S entries for this vertex
        for (int eqIdx = 0; eqIdx < numEq-dim; ++eqIdx) {
            outStream << this->curSol().base()[dofIdxGlobal*(numEq-dim) + eqIdx][0]<<" ";
        }
        // get ux, uy, uz entries for this vertex
        for (int j = 0; j< dim; ++j)
            outStream << this->curSol().base()[numScv*(numEq-dim) + dofIdxGlobal*dim + j][0] <<" ";

        int vIdxGlobal = this->dofMapper().index(entity);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vIdxGlobal);
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param entity The Entity
     *
     * Due to the mixed discretization schemes which are combined via pdelab for this model
     * the solution vector has a different form than in the pure box models
     * it sorts the primary variables in the following way:
     * p_vertex0 S_vertex0 p_vertex1 S_vertex1 p_vertex2 ....p_vertexN S_vertexN
     * ux_vertex0 uy_vertex0 uz_vertex0 ux_vertex1 uy_vertex1 uz_vertex1 ...
     *
     * Therefore, the deserializeEntity function has to be modified.
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!inStream.good()){
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << dofIdxGlobal);
        }
        int numScv = this->gridView().size(dim);
        for (int eqIdx = 0; eqIdx < numEq-dim; ++eqIdx) {
        // read p and S entries for this vertex
        inStream >> this->curSol().base()[dofIdxGlobal*(numEq-dim) + eqIdx][0];}
        for (int j = 0; j< dim; ++j){
            // read ux, uy, uz entries for this vertex
            inStream >> this->curSol().base()[numScv*(numEq-dim) + dofIdxGlobal*dim + j][0];}
    }


    /*!
     * \brief \copybrief ImplicitModel::addOutputVtkFields
     *
     * Specialization for the ElOnePTwoCBoxModel, add one-phase two-component
     * properties, solid displacement, stresses, effective properties and the
     * process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer) {
        // check whether compressive stresses are defined to be positive
        // (rockMechanicsSignConvention_ == true) or negative
        rockMechanicsSignConvention_ =  GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, RockMechanicsSignConvention);

        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // create the required scalar and vector fields
        unsigned numVertices = this->gridView_().size(dim);
        unsigned numElements = this->gridView_().size(0);

//         // create the required fields for vertex data
        ScalarField &pw = *writer.allocateManagedBuffer(numVertices);
        ScalarField &pn = *writer.allocateManagedBuffer(numVertices);
        ScalarField &pc = *writer.allocateManagedBuffer(numVertices);
        ScalarField &sw = *writer.allocateManagedBuffer(numVertices);
        ScalarField &sn = *writer.allocateManagedBuffer(numVertices);
        VectorField &displacement = *writer.template allocateManagedBuffer<Scalar, dim>(numVertices);
        ScalarField &rhoW = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rhoN = *writer.allocateManagedBuffer(numVertices);

        ScalarField &E = *writer.allocateManagedBuffer(numElements);
        ScalarField &B = *writer.allocateManagedBuffer(numElements);
        ScalarField &nu = *writer.allocateManagedBuffer(numElements);


        // create the required fields for element data
        // effective stresses
        VectorField &deltaEffStressX = *writer.template allocateManagedBuffer<Scalar,
                dim>(numElements);
        VectorField &deltaEffStressY = *writer.template allocateManagedBuffer<Scalar,
                dim>(numElements);
        VectorField &deltaEffStressZ = *writer.template allocateManagedBuffer<Scalar,
                dim>(numElements);
        // total stresses
        VectorField &totalEffStressX = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        VectorField &totalEffStressY = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        VectorField &totalEffStressZ = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        // initial stresses
        VectorField &initStressX = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        VectorField &initStressY = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        VectorField &initStressZ = *writer.template allocateManagedBuffer<
                Scalar, dim>(numElements);
        // principal stresses
        ScalarField &principalStress1 = *writer.allocateManagedBuffer(
                numElements);
        ScalarField &principalStress2 = *writer.allocateManagedBuffer(
                numElements);
        ScalarField &principalStress3 = *writer.allocateManagedBuffer(
                numElements);

        ScalarField &effKx = *writer.allocateManagedBuffer(numElements);
        ScalarField &effPorosity = *writer.allocateManagedBuffer(numElements);
        ScalarField &effPorosityOldIteration = *writer.allocateManagedBuffer(numElements);
        ScalarField &effectivePressure = *writer.allocateManagedBuffer(numElements);
        ScalarField &initialPressure = *writer.allocateManagedBuffer(numElements);
        ScalarField &deltaEffPressure = *writer.allocateManagedBuffer(numElements);

        ScalarField &Pcrtens = *writer.allocateManagedBuffer(numElements);
        ScalarField &Pcrshe = *writer.allocateManagedBuffer(numElements);

        ScalarField &Psc = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &taum = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &sigmam = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &sigmamI = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &waterContent = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &suction = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &suctionStress = *writer.allocateManagedBuffer(numElements);//khodam
        ScalarField &safetyFactor = *writer.allocateManagedBuffer(numElements);//khodam



        ScalarField &volumetricStrain = *writer.allocateManagedBuffer(numElements);

        // initialize cell stresses, cell-wise hydraulic parameters and cell pressure with zero


        for (unsigned int eIdx = 0; eIdx < numElements; ++eIdx) {
            deltaEffStressX[eIdx] = Scalar(0.0);
            if (dim >= 2)
                deltaEffStressY[eIdx] = Scalar(0.0);
            if (dim >= 3)
                deltaEffStressZ[eIdx] = Scalar(0.0);

            totalEffStressX[eIdx] = Scalar(0.0);
            if (dim >= 2)
                totalEffStressY[eIdx] = Scalar(0.0);
            if (dim >= 3)
                totalEffStressZ[eIdx] = Scalar(0.0);

            initStressX[eIdx] = Scalar(0.0);
            if (dim >= 2)
                initStressY[eIdx] = Scalar(0.0);
            if (dim >= 3)
                initStressZ[eIdx] = Scalar(0.0);

            principalStress1[eIdx] = Scalar(0.0);
            if (dim >= 2)
                principalStress2[eIdx] = Scalar(0.0);
            if (dim >= 3)
                principalStress3[eIdx] = Scalar(0.0);

            effPorosity[eIdx] = Scalar(0.0);
            effPorosityOldIteration[eIdx] = Scalar(0.0);
            effKx[eIdx] = Scalar(0.0);
            effectivePressure[eIdx] = Scalar(0.0);
            deltaEffPressure[eIdx] = Scalar(0.0);
            initialPressure[eIdx] = Scalar(0.0);

            Pcrtens[eIdx] = Scalar(0.0);
            Pcrshe[eIdx] = Scalar(0.0);

            volumetricStrain[eIdx] = Scalar(0.0);

            taum[eIdx] = Scalar(0.0);//khodam
            sigmam[eIdx] = Scalar(0.0);//khodam
            sigmamI[eIdx] = Scalar(0.0);//khodam
            waterContent[eIdx] = Scalar(0.0);//khodam
            suction[eIdx] = Scalar(0.0);//khodam
            suctionStress[eIdx] = Scalar(0.0);//khodamafetyFactor[eIdx] = Scalar(0.0);//khodam
            safetyFactor[eIdx] = Scalar(0.0);//khodam

        }

        ScalarField &rank = *writer.allocateManagedBuffer(numElements);


        FVElementGeometry fvGeometry;
        ElementVolumeVariables elemVolVars;

        const GridFunctionSpace& gridFunctionSpace = this->problem_().model().jacobianAssembler().gridFunctionSpace();
        const typename GridFunctionSpace::Ordering& ordering = gridFunctionSpace.ordering();
        // initialize start and end of element iterator
        // loop over all elements (cells)
        for (const auto& element : elements(this->gridView_())) {
            if(element.partitionType() == Dune::InteriorEntity)
            {

            // get FE function spaces to calculate gradients (gradient data of momentum balance
            // equation is not stored in fluxvars since it is not evaluated at box integration point)
            // copy the values of the sol vector to the localFunctionSpace values of the current element
            LocalFunctionSpace localFunctionSpace(gridFunctionSpace);
            localFunctionSpace.bind(element);
            std::vector<Scalar> values(localFunctionSpace.size());
            for (typename LocalFunctionSpace::Traits::IndexContainer::size_type k=0; k<localFunctionSpace.size(); ++k)
            {
                const typename GridFunctionSpace::Ordering::Traits::DOFIndex& di = localFunctionSpace.dofIndex(k);
                typename GridFunctionSpace::Ordering::Traits::ContainerIndex ci;
                ordering.mapIndex(di.view(),ci);
                values[k] = sol[ci];
            }

            // local function space for solid displacement
            const unsigned int dispSize = localFunctionSpace.child(0).size();
            typedef typename LocalFunctionSpace::template Child<0>::Type ScalarDispLFS;
            // further types required for gradient calculations
            typedef typename ScalarDispLFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
            typedef typename ScalarDispLFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;

            unsigned int eIdx = this->problem_().model().elementMapper().index(element);
            rank[eIdx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), element);
            elemVolVars.update(this->problem_(), element, fvGeometry, false);

            // loop over all local vertices of the cell
            int numScv = element.subEntities(dim);

            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
            {
                unsigned int vIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dim);

                #if PROBLEM_IS_CC==1
                pw[vIdxGlobal] = this->problem_().getpw(element, fvGeometry, 0);
                pn[vIdxGlobal] = this->problem_().getpn(element, fvGeometry, 0);
                pc[vIdxGlobal] = this->problem_().getpc(element, fvGeometry, 0);
                sw[vIdxGlobal] = this->problem_().getSw(element, fvGeometry, 0);
                sn[vIdxGlobal] = this->problem_().getSn(element, fvGeometry, 0);
                rhoW[vIdxGlobal] = this->problem_().getRhow(element, fvGeometry, 0);
                rhoN[vIdxGlobal] = this->problem_().getRhon(element, fvGeometry, 0);

                #else
                pw[vIdxGlobal] = this->problem_().getpw(element, fvGeometry, scvIdx);
                pn[vIdxGlobal] = this->problem_().getpn(element, fvGeometry, scvIdx);
                pc[vIdxGlobal] = this->problem_().getpc(element, fvGeometry, scvIdx);
                sw[vIdxGlobal] = this->problem_().getSw(element, fvGeometry, scvIdx);
                sn[vIdxGlobal] = this->problem_().getSn(element, fvGeometry, scvIdx);
                rhoW[vIdxGlobal] = this->problem_().getRhow(element, fvGeometry, scvIdx);
                rhoN[vIdxGlobal] = this->problem_().getRhon(element, fvGeometry, scvIdx);

                #endif


                // the following lines are correct for rock mechanics sign convention
                // but lead to a very counter-intuitive output therefore, they are commented.
                // in case of rock mechanics sign convention solid displacement is
                // defined to be negative if it points in positive coordinate direction
//                if(rockMechanicsSignConvention_){
//                    DimVector tmpDispl;
//                    tmpDispl = Scalar(0);
//                    tmpDispl -= elemVolVars[scvIdx].displacement();
//                    displacement[vIdxGlobal] = tmpDispl;
//                    }
//
//                else
                    displacement[vIdxGlobal] = elemVolVars[scvIdx].displacement();

                double Keff;
//                     double exponent;
//                     exponent = 22.2 * (elemVolVars[scvIdx].effPorosity
//                             / elemVolVars[scvIdx].porosity() - 1);
//                     Keff = this->problem_().spatialParams().intrinsicPermeability(*eIt, fvGeometry, scvIdx)[0][0];
//                     Keff *= exp(exponent);

                double factor;
                factor = pow( (elemVolVars[scvIdx].effPorosity / elemVolVars[scvIdx].porosity()),15);
                Keff = this->problem_().spatialParams().intrinsicPermeability(element, fvGeometry, scvIdx)[0][0];
//                     Keff *= factor;

                effKx[eIdx] += Keff;

                effectivePressure[eIdx] += (pn[vIdxGlobal] * sn[vIdxGlobal]
                                            + pw[vIdxGlobal] * sw[vIdxGlobal]);
                waterContent[eIdx] += (0.46 * sw[vIdxGlobal])/ numScv;//khodam (phi-swr)*sw+swr
                suction[eIdx] += (pn[vIdxGlobal] -  pw[vIdxGlobal])/ numScv;//khodam
                suctionStress[eIdx] += -sw[vIdxGlobal] * (pn[vIdxGlobal] -  pw[vIdxGlobal])/ numScv;//khodam

                effPorosity[eIdx] += elemVolVars[scvIdx].effPorosity;

                volumetricStrain[eIdx] += elemVolVars[scvIdx].volumetricStrain;

            }

            effKx[eIdx] = effKx[eIdx]/ numScv;

            effectivePressure[eIdx] = effectivePressure[eIdx]/ numScv;

            effPorosity[eIdx] = effPorosity[eIdx] / numScv;

            volumetricStrain[eIdx] = volumetricStrain[eIdx] / numScv;

            const auto geometry = element.geometry();

            const GlobalPosition& cellCenter = geometry.center();
            const GlobalPosition& cellCenterLocal = geometry.local(cellCenter);

            deltaEffPressure[eIdx] = effectivePressure[eIdx] + this->problem().pInit(cellCenter, cellCenterLocal, element);
            initialPressure[eIdx] = -1.0 * this->problem().pInit(cellCenter, cellCenterLocal, element);
//             pInit[eIdx] = this->problem().pInit(cellCenter, cellCenterLocal, element);
            // determin changes in effective stress from current solution
            // evaluate gradient of displacement shape functions
            std::vector<JacobianType_V> vRefShapeGradient(dispSize);
            localFunctionSpace.child(0).finiteElement().localBasis().evaluateJacobian(cellCenterLocal, vRefShapeGradient);

            // get jacobian to transform the gradient to physical element
            const JacobianInverseTransposed jacInvT = geometry.jacobianInverseTransposed(cellCenterLocal);
            std::vector < Dune::FieldVector<RF, dim> > vShapeGradient(dispSize);
            for (size_t i = 0; i < dispSize; i++) {
                vShapeGradient[i] = 0.0;
                jacInvT.umv(vRefShapeGradient[i][0], vShapeGradient[i]);
            }
            // calculate gradient of current displacement
            typedef Dune::FieldMatrix<RF, dim, dim> DimMatrix;
            DimMatrix uGradient(0.0);
            for (int coordDir = 0; coordDir < dim; ++coordDir) {
                const ScalarDispLFS & scalarDispLFS = localFunctionSpace.child(coordDir);

                for (size_t i = 0; i < scalarDispLFS.size(); i++)
                    uGradient[coordDir].axpy(values[scalarDispLFS.localIndex(i)],vShapeGradient[i]);
            }

            const Dune::FieldVector<Scalar, 3> lameParams =    this->problem_().spatialParams().lameParams(element,fvGeometry, 0);
            const Scalar lambda = lameParams[0];
            const Scalar mu = lameParams[1];

            B[eIdx] = lambda + 2.0/3.0*mu;
            E[eIdx] = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
            nu[eIdx] = lambda / (2 * (lambda + mu));

            // calculate strain tensor
            Dune::FieldMatrix<RF, dim, dim> epsilon;
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j)
                    epsilon[i][j] = 0.5 * (uGradient[i][j] + uGradient[j][i]);

            RF traceEpsilon = 0;
            for (int i = 0; i < dim; ++i)
                traceEpsilon += epsilon[i][i];

            // calculate effective stress tensor
            Dune::FieldMatrix<RF, dim, dim> sigma(0.0);
            for (int i = 0; i < dim; ++i) {
                sigma[i][i] = lambda * traceEpsilon;
                for (int j = 0; j < dim; ++j)
                    sigma[i][j] += 2.0 * mu * epsilon[i][j];
            }

            // in case of rock mechanics sign convention compressive stresses
            // are defined to be positive
            if(rockMechanicsSignConvention_){
                deltaEffStressX[eIdx] -= sigma[0];
                if (dim >= 2) {
                    deltaEffStressY[eIdx] -= sigma[1];
                }
                if (dim >= 3) {
                    deltaEffStressZ[eIdx] -= sigma[2];
                }
            }
            else{
                deltaEffStressX[eIdx] = sigma[0];
                if (dim >= 2) {
                    deltaEffStressY[eIdx] = sigma[1];
                }
                if (dim >= 3) {
                    deltaEffStressZ[eIdx] = sigma[2];
                }
            }

            // retrieve prescribed initial stresses from problem file
            DimVector tmpInitStress = this->problem_().initialStress(cellCenter, 0);
            if(rockMechanicsSignConvention_){
                initStressX[eIdx][0] = tmpInitStress[0];
                if (dim >= 2) {
                    initStressY[eIdx][1] = tmpInitStress[1];
                    }
                if (dim >= 3) {
                    initStressZ[eIdx][2] = tmpInitStress[2];
                }
            }
            else{
                initStressX[eIdx][0] -= tmpInitStress[0];
                if (dim >= 2) {
                    initStressY[eIdx][1] -= tmpInitStress[1];
                    }
                if (dim >= 3) {
                    initStressZ[eIdx][2] -= tmpInitStress[2];
                }
            }

            // calculate total stresses
            // in case of rock mechanics sign convention compressive stresses
            // are defined to be positive and total stress is calculated by adding the pore pressure
            if(rockMechanicsSignConvention_){
                totalEffStressX[eIdx][0] = initStressX[eIdx][0] + deltaEffStressX[eIdx][0]    - initialPressure[eIdx];
                if (dim >= 2) {
                    totalEffStressX[eIdx][1] = initStressX[eIdx][1] + deltaEffStressX[eIdx][1];
                    totalEffStressY[eIdx][0] = initStressY[eIdx][0] + deltaEffStressY[eIdx][0];
                    totalEffStressY[eIdx][1] = initStressY[eIdx][1] + deltaEffStressY[eIdx][1]    - initialPressure[eIdx];
                }
                if (dim >= 3) {
                    totalEffStressX[eIdx][2] = initStressX[eIdx][2] + deltaEffStressX[eIdx][2];
                    totalEffStressY[eIdx][2] = initStressY[eIdx][2] + deltaEffStressY[eIdx][2];
                    totalEffStressZ[eIdx][0] = initStressZ[eIdx][0] + deltaEffStressZ[eIdx][0];
                    totalEffStressZ[eIdx][1] = initStressZ[eIdx][1] + deltaEffStressZ[eIdx][1];
                    totalEffStressZ[eIdx][2] = initStressZ[eIdx][2] + deltaEffStressZ[eIdx][2]    - initialPressure[eIdx];
                }
            }
            else{
                totalEffStressX[eIdx][0] = initStressX[eIdx][0] + deltaEffStressX[eIdx][0]    + initialPressure[eIdx];
                if (dim >= 2) {
                    totalEffStressX[eIdx][1] = initStressX[eIdx][1] + deltaEffStressX[eIdx][1];
                    totalEffStressY[eIdx][0] = initStressY[eIdx][0] + deltaEffStressY[eIdx][0];
                    totalEffStressY[eIdx][1] = initStressY[eIdx][1] + deltaEffStressY[eIdx][1]    + initialPressure[eIdx];
                }
                if (dim >= 3) {
                    totalEffStressX[eIdx][2] = initStressX[eIdx][2] + deltaEffStressX[eIdx][2];
                    totalEffStressY[eIdx][2] = initStressY[eIdx][2] + deltaEffStressY[eIdx][2];
                    totalEffStressZ[eIdx][0] = initStressZ[eIdx][0] + deltaEffStressZ[eIdx][0];
                    totalEffStressZ[eIdx][1] = initStressZ[eIdx][1] + deltaEffStressZ[eIdx][1];
                    totalEffStressZ[eIdx][2] = initStressZ[eIdx][2] + deltaEffStressZ[eIdx][2]    + initialPressure[eIdx];
                }
            }
        }
        }

        // calculate principal stresses i.e. the eigenvalues of the total stress tensor
        Scalar a1, a2, a3;
        DimMatrix totalEffStress;
        DimVector eigenValues;

        for (unsigned int eIdx = 0; eIdx < numElements; eIdx++)
        {
            eigenValues = Scalar(0);
            totalEffStress = Scalar(0);

            totalEffStress[0] = totalEffStressX[eIdx];
            if (dim >= 2)
                totalEffStress[1] = totalEffStressY[eIdx];
            if (dim >= 3)
                totalEffStress[2] = totalEffStressZ[eIdx];

            calculateEigenValues<dim>(eigenValues, totalEffStress);


            for (int i = 0; i < dim; i++)
                {
                    if (std::isnan(eigenValues[i]))
                        eigenValues[i] = 0.0;
                }

            // sort principal stresses: principalStress1 >= principalStress2 >= principalStress3
            if (dim == 2) {
                a1 = eigenValues[0];
                a2 = eigenValues[1];


                if (a1 >= a2) {
                    principalStress1[eIdx] = a1;
                    principalStress2[eIdx] = a2;
                    principalStress3[eIdx] = a2;//khodam for dim == 2
                } else {
                    principalStress1[eIdx] = a2;
                    principalStress2[eIdx] = a1;
                    principalStress3[eIdx] = a1;//khodam for dim == 2
                }
            }

            if (dim == 3) {
                a1 = eigenValues[0];
                a2 = eigenValues[1];
                a3 = eigenValues[2];

                if (a1 >= a2) {
                    if (a1 >= a3) {
                        principalStress1[eIdx] = a1;
                        if (a2 >= a3) {
                            principalStress2[eIdx] = a2;
                            principalStress3[eIdx] = a3;
                        }
                        else //a3 > a2
                        {
                            principalStress2[eIdx] = a3;
                            principalStress3[eIdx] = a2;
                        }
                    }
                    else // a3 > a1
                    {
                        principalStress1[eIdx] = a3;
                        principalStress2[eIdx] = a1;
                        principalStress3[eIdx] = a2;
                    }
                } else // a2>a1
                {
                    if (a2 >= a3) {
                        principalStress1[eIdx] = a2;
                        if (a1 >= a3) {
                            principalStress2[eIdx] = a1;
                            principalStress3[eIdx] = a3;
                        }
                        else //a3>a1
                        {
                            principalStress2[eIdx] = a3;
                            principalStress3[eIdx] = a1;
                        }
                    }
                    else //a3>a2
                    {
                        principalStress1[eIdx] = a3;
                        principalStress2[eIdx] = a2;
                        principalStress3[eIdx] = a1;
                    }
                }
            }
//             Scalar taum  = 0.0;
//             Scalar sigmam = 0.0;
            Scalar Peff = effectivePressure[eIdx];



for (const auto& element : elements(this->gridView_())) {//khodam
if(element.partitionType() == Dune::InteriorEntity)//khodam
{
int numScv = element.subEntities(dim);//khodam

for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)//khodam
    {
            unsigned int vIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dim);

//             Scalar theta = M_PI / 6;
//             Scalar S0 = 0.0;
            Scalar theta = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.FrictionAngle)*M_PI / 180.0 ; //khodam
            Scalar S0 = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.Cohesion);//khodam

            using std::abs;
            using std::sin;
            using std::cos;

            taum[eIdx] = (principalStress1[eIdx] - principalStress3[eIdx]) / 2; //sigmaII
            sigmam[eIdx] = (principalStress1[eIdx] + principalStress3[eIdx]) / 2;
            sigmamI[eIdx] = abs(sigmam[eIdx]) - suctionStress[eIdx];//khodam
//             sigmamI[eIdx] = (abs(principalStress3[eIdx] + principalStress1[eIdx]) / 2.0) + sw[vIdxGlobal] * (pn[vIdxGlobal] - pw[vIdxGlobal]);//khodam after Greco paper

            Scalar Psc = -abs(taum[eIdx]) / sin(theta) + S0 * cos(theta) / sin(theta)
                    + sigmam[eIdx]; //critical fluid pressure for shear failure
            // Pressure margins according to J. Rutqvist et al. / International Journal of Rock Mecahnics & Mining Sciences 45 (2008), 132-143
            Pcrtens[eIdx] = Peff - principalStress3[eIdx];
            Pcrshe[eIdx] = Peff - Psc;//equation 2.61 mellanie thesis
            //LFS khodam
            if (std::abs(taum[eIdx]) < 1.0e-2){
                safetyFactor[eIdx] = 5.0;
             }
            else
            {

//                  safetyFactor[eIdx] = cos(theta) / fabs(taum) * (S0 + sigmamI * tan(theta));//khodam
                safetyFactor[eIdx] = cos(theta) / (taum[eIdx]) * (S0 + (sigmamI[eIdx]) * tan(theta));//khodam

                 if (safetyFactor[eIdx]<1)
                      safetyFactor[eIdx]= 1;
                 else if (safetyFactor[eIdx]>3)
                     safetyFactor[eIdx]= 3;
             }

        }
}
}
        }

        writer.attachVertexData(pw, "pw");
        writer.attachVertexData(pn, "pN");
        writer.attachVertexData(pc, "pC");
        writer.attachVertexData(sw, "SW");
        writer.attachVertexData(sn, "SN");
        writer.attachVertexData(rhoW, "rhoW");
        writer.attachVertexData(rhoN, "rhoN");
        writer.attachVertexData(displacement, "u", dim);

        writer.attachCellData(deltaEffStressX, "effectivestresschanges X", dim);
        if (dim >= 2)
            writer.attachCellData(deltaEffStressY, "effectivestresschanges Y",    dim);
        if (dim >= 3)
            writer.attachCellData(deltaEffStressZ, "effectivestresschanges Z",    dim);

        writer.attachCellData(principalStress1, "principalstress1");
        if (dim >= 2)
            writer.attachCellData(principalStress2, "principalstress2");
//         if (dim >= 3)
//             writer.attachCellData(principalStress3, "principalstress3");

        writer.attachCellData(totalEffStressX, "totalEffStressesX", dim);
        if (dim >= 2)
            writer.attachCellData(totalEffStressY, "totalEffStressesY", dim);
        if (dim >= 3)
            writer.attachCellData(totalEffStressZ, "totalEffStressesZ", dim);

        writer.attachCellData(initStressX, "initialstressesX", dim);
        if (dim >= 2)
            writer.attachCellData(initStressY, "initialstressesY", dim);
        if (dim >= 3)
            writer.attachCellData(initStressZ, "initialstressesZ", dim);

        writer.attachCellData(deltaEffPressure, "deltapEff");
        writer.attachCellData(effectivePressure, "effectivePressure");
        writer.attachCellData(waterContent, "waterContent");//khodam
        writer.attachCellData(suction, "suction");//khodam
        writer.attachCellData(suctionStress, "suctionStress");//khodam
        writer.attachCellData(initialPressure, "pInit");
        writer.attachCellData(Pcrtens, "Pcr_tensile");
        writer.attachCellData(Pcrshe, "Pcrshe");
        writer.attachCellData(effKx, "effectiveKxx");
        writer.attachCellData(effPorosity, "effectivePorosity");

        writer.attachCellData(E, "E");
        writer.attachCellData(B, "B");
        writer.attachCellData(nu, "nu");
        writer.attachCellData(Psc, "Psc");//khodam
        writer.attachCellData(taum, "taum");//khodam
        writer.attachCellData(sigmam, "sigmam");//khodam
        writer.attachCellData(sigmamI, "sigmamI");//khodam

        writer.attachCellData(volumetricStrain, "volumetricStrain");
        writer.attachCellData(safetyFactor, "FactorOfSafety");//khodam
    }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     */
    void applyInitialSolution_() {
//         typedef typename GET_PROP_TYPE(TypeTag, InitialPressSat) InitialPressSat;
//         InitialPressSat initialPressSat(this->problem_().gridView());
//         std::cout << "el2pmodel calls: initialPressSat" << std::endl;
//         initialPressSat.setPressure(this->problem_().pInit());
//
//         typedef typename GET_PROP_TYPE(TypeTag, InitialDisplacement) InitialDisplacement;
//         InitialDisplacement initialDisplacement(this->problem_().gridView());
//
//         typedef Dune::PDELab::CompositeGridFunction<InitialPressSat,
//                 InitialDisplacement> InitialSolution;
//         InitialSolution initialSolution(initialPressSat, initialDisplacement);
//
//         int numDofs = this->jacobianAssembler().gridFunctionSpace().size();
//         //this->curSol().resize(numDofs);
//         //this->prevSol().resize(numDofs);
//         std::cout << "numDofs = " << numDofs << std::endl;
//
//         Dune::PDELab::interpolate(initialSolution,
//                 this->jacobianAssembler().gridFunctionSpace(), this->curSol());
//         Dune::PDELab::interpolate(initialSolution,
//                 this->jacobianAssembler().gridFunctionSpace(), this->prevSol());
    }

    const Problem& problem() const {
        return this->problem_();
    }

private:
    bool rockMechanicsSignConvention_;
    Scalar pnRef_ = 1e5;

};
}
#include "propertydefaults.hh"
#endif
