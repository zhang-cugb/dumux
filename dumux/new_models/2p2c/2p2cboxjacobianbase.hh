/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH
#define DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/2p2c/2p2ctraits.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

//// forward declaration of the 2p2c box model
//template<class ProblemT, class TwoPTwoCTraitsT>
//class TwoPTwoCBoxModel;

///////////////////////////////////////////////////////////////////////////
// TwoPTwoCBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief 2P-2C specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxJacobian for the 2P-2C twophase flow.
 */
template<class ProblemT,
         class BoxTraitsT,
         class TwoPTwoCTraitsT,
         class Implementation>
class TwoPTwoCBoxJacobianBase : public BoxJacobian<ProblemT,
                                                   BoxTraitsT,
                                                   Implementation,
                                                   typename TwoPTwoCTraitsT::ElementData,
                                                   typename TwoPTwoCTraitsT::VertexData>

{
protected:
    typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                    BoxTraitsT,
                                    TwoPTwoCTraitsT,
                                    Implementation>              ThisType;
    typedef BoxJacobian<ProblemT, 
                        BoxTraitsT,
                        Implementation,
                        typename TwoPTwoCTraitsT::ElementData,
                        typename TwoPTwoCTraitsT::VertexData>    ParentType;
    
    typedef ProblemT                                Problem;
    typedef typename Problem::DomainTraits          DomTraits;
    typedef BoxTraitsT                              BoxTraits;
    typedef TwoPTwoCTraitsT                         TwoPTwoCTraits;

    enum {
        dim              = DomTraits::dim,
        dimWorld         = DomTraits::dimWorld,

        numEq            = BoxTraits::numEq,
        numPhases        = TwoPTwoCTraits::numPhases,
        numComponents    = TwoPTwoCTraits::numComponents,

        pressureIdx      = TwoPTwoCTraits::pressureIdx,
        switchIdx        = TwoPTwoCTraits::switchIdx,

        wPhase           = TwoPTwoCTraits::wPhase,
        nPhase           = TwoPTwoCTraits::nPhase,

        wComp            = TwoPTwoCTraits::wComp,
        nComp            = TwoPTwoCTraits::nComp,

        wPhaseOnly       = TwoPTwoCTraits::wPhaseOnly,
        nPhaseOnly       = TwoPTwoCTraits::nPhaseOnly,
        bothPhases       = TwoPTwoCTraits::bothPhases
    };
    static const int formulation  = TwoPTwoCTraits::formulation;
    enum {
        pWsN             = TwoPTwoCTraits::pWsN,
        pNsW             = TwoPTwoCTraits::pNsW,
    };


    typedef typename DomTraits::Scalar                Scalar;
    typedef typename DomTraits::CoordScalar           CoordScalar;
    typedef typename DomTraits::Grid                  Grid;
    typedef typename DomTraits::Vertex                Vertex;
    typedef typename DomTraits::Element               Element;
    typedef typename DomTraits::ElementIterator       ElementIterator;
    typedef typename Element::EntityPointer           ElementPointer;
    typedef typename DomTraits::LocalPosition         LocalPosition;
    typedef typename DomTraits::GlobalPosition        GlobalPosition;
    typedef typename DomTraits::VertexIterator        VertexIterator;

    typedef typename BoxTraits::SolutionVector      SolutionVector;
    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;
    typedef typename Grid::CollectiveCommunication  CollectiveCommunication;

    typedef typename TwoPTwoCTraits::PhasesVector        PhasesVector;
    typedef typename TwoPTwoCTraits::ElementData         ElementData;
    typedef typename TwoPTwoCTraits::VertexData          VertexData;

    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVertexData {
        int phaseState;
        int oldPhaseState;
    };

public:
    TwoPTwoCBoxJacobianBase(ProblemT &problem)
        : ParentType(problem),
          staticVertexDat_(problem.numVertices())
    {
        switchFlag_ = false;
    };

    /*!
     * \brief Function to update variable data of the vertices of the
     *        the current element (essentially all secondary variables)
     */
    void updateVertexData_(ElementData &elemDat,
                           const LocalFunction &sol,
                           int vertIdx,
                           bool isOldSol) const
    {
        const Element       &element = this->curElement_();
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            DomTraits::referenceElement(element.type()).position(vertIdx,
                                                                 dim);
        
        VertexData           &vertDat = elemDat[vertIdx];
        const SolutionVector &vertSol = sol[vertIdx];
        int globalVertIdx = this->problem_.vertexIdx(this->curElement_(), 
                                                     vertIdx);

        Scalar temperature = asImp_()->temperature(vertSol);

        int phaseState;
        if (isOldSol)
            phaseState = staticVertexDat_[globalVertIdx].oldPhaseState;
        else
            phaseState = staticVertexDat_[globalVertIdx].phaseState;

        if (formulation == pWsN)
        {
            vertDat.pressure[wPhase] = vertSol[pressureIdx];
            if (phaseState == bothPhases) vertDat.saturation[nPhase] = vertSol[switchIdx];
            else if (phaseState == wPhaseOnly) vertDat.saturation[nPhase] = 0.0;
            else if (phaseState == nPhaseOnly) vertDat.saturation[nPhase] = 1.0;
            else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

            vertDat.saturation[wPhase] = 1.0 - vertDat.saturation[nPhase];
            vertDat.pC = this->problem_.materialLaw().pC(vertDat.saturation[wPhase],
                                                         global,
                                                         element,
                                                         local);
            
            vertDat.pressure[nPhase] = vertDat.pressure[wPhase] + vertDat.pC;
        }
        else if (formulation == pNsW)
        {
            vertDat.pressure[nPhase] = vertSol[pressureIdx];
            if (phaseState == bothPhases) vertDat.saturation[wPhase] = vertSol[switchIdx];
            else if (phaseState == wPhaseOnly) vertDat.saturation[wPhase] = 1.0;
            else if (phaseState == nPhaseOnly) vertDat.saturation[wPhase] = 0.0;
            else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

            vertDat.saturation[nPhase] = 1.0 - vertDat.saturation[wPhase];
            vertDat.pC = this->problem_.materialLaw().pC(vertDat.saturation[wPhase],
                                                         global,
                                                         element,
                                                         local);

            vertDat.pressure[wPhase] = vertDat.pressure[nPhase] - vertDat.pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        // Solubilities of components in phases
        if (phaseState == bothPhases) {
            vertDat.massfrac[wComp][nPhase] = this->problem_.multicomp().xWN(vertDat.pressure[nPhase], temperature);
            vertDat.massfrac[nComp][wPhase] = this->problem_.multicomp().xAW(vertDat.pressure[nPhase], temperature);
        }
        else if (phaseState == wPhaseOnly) {
            vertDat.massfrac[wComp][nPhase] = 0.0;
            vertDat.massfrac[nComp][wPhase] = vertSol[switchIdx];
        }
        else if (phaseState == nPhaseOnly){
            vertDat.massfrac[wComp][nPhase] = vertSol[switchIdx];
            vertDat.massfrac[nComp][wPhase] = 0.0;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

        vertDat.massfrac[wComp][wPhase] = 1.0 - vertDat.massfrac[nComp][wPhase];
        vertDat.massfrac[nComp][nPhase] = 1.0 - vertDat.massfrac[wComp][nPhase];

        // Densities
        vertDat.density[wPhase] = this->problem_.wettingPhase().density(temperature,
                                                                        vertDat.pressure[wPhase],
                                                                        vertDat.massfrac[nComp][wPhase]);
        vertDat.density[nPhase] = this->problem_.nonwettingPhase().density(temperature,
                                                                           vertDat.pressure[nPhase],
                                                                           vertDat.massfrac[wComp][nPhase]);

        // Mobilities
        vertDat.mobility[wPhase] = this->problem_.materialLaw().mobW(vertDat.saturation[wPhase],
                                                                     global,
                                                                     element,
                                                                     local,
                                                                     temperature,
                                                                     vertDat.pressure[wPhase]);
        vertDat.mobility[nPhase] = this->problem_.materialLaw().mobN(vertDat.saturation[nPhase],
                                                                     global,
                                                                     element,
                                                                     local,
                                                                     temperature,
                                                                     vertDat.pressure[nPhase]);

        // diffusion coefficents
        vertDat.diffCoeff[wPhase] = this->problem_.wettingPhase().diffCoeff(temperature, vertDat.pressure[wPhase]);
        vertDat.diffCoeff[nPhase] = this->problem_.nonwettingPhase().diffCoeff(temperature, vertDat.pressure[nPhase]);

        // porosity
        vertDat.porosity = this->problem_.soil().porosity(global,
                                                          this->curElement_(), 
                                                          local);

    }

    struct FluxVars
    {
        typedef typename FVElementGeometry::SubControlVolume     SCV;
        typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

        FluxVars(const FVElementGeometry &elemGeom,
                 int faceIdx)
            : fvElemGeom(elemGeom)
        {
            face = &fvElemGeom.subContVolFace[faceIdx];

            int i = face->i;
            int j = face->j;
            
            insideSCV = &fvElemGeom.subContVol[i];
            outsideSCV = &fvElemGeom.subContVol[j];

            temperatureGrad = Scalar(0.0);
            for (int phase = 0; phase < numPhases; ++phase) {
                pressureGrad[phase] = Scalar(0);
                concentrationGrad[phase] = Scalar(0);
            }
        };
        
        void calculateGradients(const Problem &problem,
                                const Element &element,
                                const ElementData &elemDat)
        {
            // calculate gradients
            GlobalPosition tmp(0.0);
            for (int idx = 0; 
                 idx < fvElemGeom.numVertices;
                 idx++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const LocalPosition &feGrad = face->grad[idx];

                // compute sum of pressure gradients for each phase
                for (int phase = 0; phase < numPhases; phase++)
                {
                    // the pressure gradient
                    tmp = feGrad;
                    tmp *= elemDat[idx].pressure[phase];
                    pressureGrad[phase] += tmp;

                    // phase density
                    densityAtIP[phase] 
                        += 
                        elemDat[idx].density[phase] *
                        face->shapeValue[idx];
                }

                // the concentration gradient of the non-wetting
                // component in the wetting phase
                tmp = feGrad;
                tmp *= elemDat[idx].massfrac[nComp][wPhase];
                concentrationGrad[wPhase] += tmp;

                // the concentration gradient of the wetting component
                // in the non-wetting phase
                tmp = feGrad;
                tmp *= elemDat[idx].massfrac[wComp][nPhase];
                concentrationGrad[nPhase] += tmp;

                // temperature gradient
                // TODO
                /*
                jac.asImp_().updateTempGrad(temperatureGrad,
                                            feGrad, 
                                            jac.curSol_, 
                                            idx);
                */
            }

            // correct the pressure gradients by the hydrostatic
            // pressure due to gravity
            for (int phase=0; phase < numPhases; phase++)
            {
                tmp = problem.gravity();
                tmp *= densityAtIP[phase];

                pressureGrad[phase] -= tmp;
            }
        }

        void calculateVelocities(const Problem &problem,
                                 const Element &element,
                                 const ElementData &elemDat)
        {
            // calculate the permeability tensor
            Tensor K         = problem.soil().K(insideSCV->global,
                                                element,
                                                insideSCV->local);
            const Tensor &Kj = problem.soil().K(outsideSCV->global, 
                                                element, 
                                                outsideSCV->local);
            harmonicMeanK_(K, Kj);
            
            // temporary vector for the Darcy velocity
            GlobalPosition vDarcy;
            for (int phase=0; phase < numPhases; phase++)
            {
                K.mv(pressureGrad[phase], vDarcy);  // vDarcy = K * grad p
                vDarcyNormal[phase] = vDarcy*face->normal;
            }

            // upstream and downstream vertices
            for (int phase = 0; phase < numPhases; ++phase)
            {
                upstreamIdx[phase] = face->i;
                downstreamIdx[phase] = face->j;
                if (vDarcyNormal[phase] > 0) {
                    std::swap(upstreamIdx[phase],
                              downstreamIdx[phase]);
                }
            }
        }

        void calculateDiffCoeffPM(const Problem &problem,
                                  const Element &element,
                                  const ElementData &elemDat)
        {
            const VertexData &vDat_i = elemDat[face->i];
            const VertexData &vDat_j = elemDat[face->j];

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                // make sure to only calculate diffusion coefficents
                // for phases which exist in both finite volumes
                if (vDat_i.saturation[phaseIdx] == 0 ||
                    vDat_j.saturation[phaseIdx] == 0) 
                {
                    diffCoeffPM[phaseIdx] = 0.0;
                    continue;
                }
                
                // calculate tortuosity at the nodes i and j needed
                // for porous media diffusion coefficient
                Scalar tau_i = 
                    1.0/vDat_i.porosity *
                    pow(vDat_i.porosity * vDat_i.saturation[phaseIdx], 7.0/3);
                Scalar tau_j = 
                    1.0/vDat_j.porosity *
                    pow(vDat_j.porosity * vDat_j.saturation[phaseIdx], 7.0/3);
        
                // Diffusion coefficient in the porous medium
                
                // -> arithmetic mean
                diffCoeffPM[phaseIdx]
                    = 1./2*(vDat_i.porosity * vDat_i.saturation[phaseIdx] * tau_i * vDat_i.diffCoeff[phaseIdx] +
                            vDat_j.porosity * vDat_j.saturation[phaseIdx] * tau_j * vDat_j.diffCoeff[phaseIdx]);
                // -> harmonic mean
                // = harmonicMean_(vDat_i.porosity * vDat_i.saturation[phaseIdx] * tau_i * vDat_i.diffCoeff[phaseIdx],
                //                 vDat_j.porosity * vDat_j.saturation[phaseIdx] * tau_j * vDat_j.diffCoeff[phaseIdx]);
            }
        }

        const FVElementGeometry &fvElemGeom;
        const SCVFace *face;
        const SCV     *insideSCV;
        const SCV     *outsideSCV;
        
        // gradients
        GlobalPosition pressureGrad[numPhases];
        GlobalPosition concentrationGrad[numPhases];
        GlobalPosition temperatureGrad;

        // density of each face at the integration point
        PhasesVector densityAtIP;
        
        // darcy velocity in direction of the face normal
        PhasesVector vDarcyNormal;

        // Upwind parameter
        static const Scalar upwindAlpha = 1.0; // -> use only the upstream vertex

        // local index of the upwind vertex for each phase
        int upstreamIdx[numPhases];
        // local index of the downwind vertex for each phase
        int downstreamIdx[numPhases];

        // the diffusion coefficient for the porous medium
        PhasesVector diffCoeffPM;
    };


    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(SolutionVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementData &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];
        
        // storage of component water
        result[wComp] =
            vertDat.porosity*(vertDat.density[wPhase]*
                              vertDat.saturation[wPhase]*
                              vertDat.massfrac[wComp][wPhase]
                              +
                              vertDat.density[nPhase]*
                              vertDat.saturation[nPhase]*
                              vertDat.massfrac[wComp][nPhase]);

        // storage of component air
        result[nComp] =
            vertDat.porosity*(vertDat.density[nPhase]*
                              vertDat.saturation[nPhase]*
                              vertDat.massfrac[nComp][nPhase]
                              + 
                              vertDat.density[wPhase]*
                              vertDat.saturation[wPhase]*
                              vertDat.massfrac[nComp][wPhase]);
        
        // storage of energy (if nonisothermal model is used)
        asImp_()->heatStorage(result, scvIdx, elemDat);
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(SolutionVector &flux, int faceIdx) const
    {
        FluxVars vars(this->curElementGeom_, faceIdx);
        vars.calculateGradients(this->problem_,
                                this->curElement_(),
                                this->curElemDat_);
        vars.calculateVelocities(this->problem_,
                                 this->curElement_(),
                                 this->curElemDat_);
        vars.calculateDiffCoeffPM(this->problem_,
                                  this->curElement_(),
                                  this->curElemDat_);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);
    }
    
    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeAdvectiveFlux(SolutionVector &flux, 
                              FluxVars       &vars) const
    {
        // required variables of the upstream and the downstream
        // vertices
        const VertexData &upW = this->curElemDat_[vars.upstreamIdx[wPhase]];
        const VertexData &dnW = this->curElemDat_[vars.downstreamIdx[wPhase]];

        const VertexData &upN = this->curElemDat_[vars.upstreamIdx[nPhase]];
        const VertexData &dnN = this->curElemDat_[vars.downstreamIdx[nPhase]];

        ////////
        // advective flux of the wetting component
        ////////

        // flux in the wetting phase
        flux[wComp] +=  
            vars.vDarcyNormal[wPhase] * (
                vars.upwindAlpha* // upstream vertex
                (  upW.density[wPhase] *
                   upW.mobility[wPhase] *
                   upW.massfrac[wComp][wPhase])
                +
                (1 - vars.upwindAlpha)* // downstream vertex
                (  dnW.density[wPhase] *
                   dnW.mobility[wPhase] *
                   dnW.massfrac[wComp][wPhase]));
        
        // flux in the non-wetting phase
        flux[wComp] += 
            vars.vDarcyNormal[nPhase] * (
                vars.upwindAlpha* // upstream vertex
                (  upN.density[nPhase] *
                   upN.mobility[nPhase] *
                   upN.massfrac[wComp][nPhase])
                +
                (1 - vars.upwindAlpha)* // downstream vertex
                (  dnN.density[nPhase] *
                   dnN.mobility[nPhase] *
                   dnN.massfrac[wComp][nPhase]) );
        
        ////////
        // advective flux of the non-wetting component
        ////////

        // flux in the non-wetting phase
        flux[nComp]  += 
            vars.vDarcyNormal[nPhase] * (
                vars.upwindAlpha * // upstream vertex
                (  upN.density[nPhase] *
                   upN.mobility[nPhase] *
                   upN.massfrac[nComp][nPhase])
                +
                (1 - vars.upwindAlpha) * // downstream vertex
                (  dnN.density[nPhase] *
                   dnN.mobility[nPhase] *
                   dnN.massfrac[nComp][nPhase]) );
        
        // flux in the wetting phase
        flux[nComp] += 
            vars.vDarcyNormal[wPhase] * (
                vars.upwindAlpha * // upstream vertex
                (  upW.density[wPhase] *
                   upW.mobility[wPhase] *
                   upW.massfrac[nComp][wPhase])
                +
                (1 - vars.upwindAlpha) * // downstream vertex
                (  dnW.density[wPhase] *
                   dnW.mobility[wPhase] *
                   dnW.massfrac[nComp][wPhase]) );
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeDiffusiveFlux(SolutionVector &flux, FluxVars &vars) const
    {        
        // add diffusive flux of non-wetting component in wetting phase
        Scalar tmp = vars.diffCoeffPM[wPhase] * vars.densityAtIP[wPhase] * 
            (vars.concentrationGrad[wPhase]*vars.face->normal);
        flux[nComp] += tmp; 
        flux[wComp] -= tmp;
        
        // add diffusive flux of wetting component in non-wetting phase
        tmp = vars.diffCoeffPM[nPhase] * vars.densityAtIP[nPhase] * 
            (vars.concentrationGrad[nPhase]*vars.face->normal);;
        flux[wComp] += tmp;
        flux[nComp] -= tmp;

        // TODO: the diffusive flux of the wetting component in the
        // wetting phase does rarly exhibit the same mass as the flux
        // of the non-wetting component, which means that it is not
        // -tmp
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(SolutionVector &q, int localVertexIdx)
    {
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }


    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * Called by TwoPTwoCBoxModel::initial()
     */
    void initStaticData()
    {
        setSwitched(false);

        VertexIterator it = this->problem_.vertexBegin();
        VertexIterator endit = this->problem_.vertexEnd();
        for (; it != endit; ++it)
            {
                int globalIdx = this->problem_.vertexIdx(*it);
                const GlobalPosition &globalPos = it->geometry().corner(0);

                // initialize phase state
                staticVertexDat_[globalIdx].phaseState =
                    this->problem_.initialPhaseState(*it, globalIdx, globalPos);
                staticVertexDat_[globalIdx].oldPhaseState =
                    staticVertexDat_[globalIdx].phaseState;
            }
    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     */
    void updateStaticData(SpatialFunction &curGlobalSol, SpatialFunction &oldGlobalSol)
    {
        bool wasSwitched = false;

        VertexIterator it = this->problem_.vertexBegin();
        for (; it != this->problem_.vertexEnd(); ++it)
            {
                int globalIdx = this->problem_.vertexIdx(*it);
                const GlobalPosition &global = it->geometry().corner(0);

                wasSwitched = primaryVarSwitch_(curGlobalSol,
                                                globalIdx,
                                                global)
                    || wasSwitched;
            }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        wasSwitched = this->problem_.grid().comm().max(wasSwitched);

        setSwitched(wasSwitched);
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhaseState()
    {
        int numVertices = this->problem_.numVertices();
        for (int i = 0; i < numVertices; ++i)
            staticVertexDat_[i].oldPhaseState = staticVertexDat_[i].phaseState;
    }

    /*!
     * \brief Reset the current phase state of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhaseState()
    {
        int numVertices = this->problem_.numVertices();
        for (int i = 0; i < numVertices; ++i)
            staticVertexDat_[i].phaseState = staticVertexDat_[i].oldPhaseState;
    }

    /*!
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Set whether there was a primary variable switch after in the last
     *        timestep.
     */
    void setSwitched(bool yesno)
    {
        switchFlag_ = yesno;
    }

    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SpatialFunction &globalSol, Dune::FieldVector<Scalar, 4> &mass)
    {
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        unsigned numVertices = this->problem_.numVertices();
        LocalFunction curSol(numVertices);
        ElementData   elemDat;
        VertexData tmp;
        int state;
        Scalar vol, poro, rhoN, rhoW, satN, satW, xAW, xWW, xWN, xAN, pW, Te;
        Scalar massNComp(0.), massNCompNPhase(0.), massWComp(0.), massWCompWPhase(0.);

        enum
        {   gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase state

        mass = 0;
        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
        Scalar minTe = 1e100;
        Scalar maxTe = -1e100;
        Scalar minX = 1e100;
        Scalar maxX = -1e100;

        // Loop over elements
        for (; elementIt != endit; ++elementIt)
        {

            setCurrentElement(*elementIt);
            this->restrictToElement(curSol, globalSol);
            updateElementData_(elemDat, curSol, false);
            // get geometry type

            int numLocalVerts = elementIt->template count<dim>();

            // Loop over element vertices
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                vol = this->curElementGeom_.subContVol[i].volume;

                state =  staticVertexDat_[globalIdx].phaseState;
                poro = this->problem_.porosity(this->curElement_(), i);
                rhoN = elemDat[i].density[nPhase];
                rhoW = elemDat[i].density[wPhase];
                satN = elemDat[i].saturation[nPhase];
                satW = elemDat[i].saturation[wPhase];
                xAW = elemDat[i].massfrac[nComp][wPhase];
                xWW = elemDat[i].massfrac[wComp][wPhase];
                xWN = elemDat[i].massfrac[wComp][nPhase];
                xAN = elemDat[i].massfrac[nComp][nPhase];
                pW = elemDat[i].pressure[wPhase];
                Te = Implementation::temperature_((*globalSol)[globalIdx]);
                massNComp = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);
                massNCompNPhase = vol * poro * satN * rhoN * xAN;
                massWComp = vol * poro * (satW * rhoW * xWW + satN * rhoN * xWN);
                massWCompWPhase = vol * poro * satW * rhoW * xWW;

                // get minimum and maximum values of primary variables
                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
                minX = std::min(minX, xAW);
                maxX = std::max(maxX, xAW);
                minTe = std::min(minTe, Te);
                maxTe = std::max(maxTe, Te);

                // calculate total mass
                mass[0] += massNComp;       // total mass of nonwetting component
                mass[1] += massNCompNPhase; // mass of nonwetting component in nonwetting phase
                mass[2] += massWComp;       // total mass of wetting component
                mass[3] += massWCompWPhase; // mass of wetting component in wetting phase
            }
        }
        
        // IF PARALLEL: calculate total mass including all processors
        // also works for sequential calculation
        mass = this->problem_.grid().comm().sum(mass);
        
        if(this->problem_.grid().comm() == 0) // IF PARALLEL: only print by processor with rank() == 0
        {
            // print minimum and maximum values
            std::cout << "nonwetting phase saturation: min = "<< minSat
                      << ", max = "<< maxSat << std::endl;
            std::cout << "wetting phase pressure: min = "<< minP
                      << ", max = "<< maxP << std::endl;
            std::cout << "mass fraction nComp: min = "<< minX
                      << ", max = "<< maxX << std::endl;
            std::cout << "temperature: min = "<< minTe
                      << ", max = "<< maxTe << std::endl;
        }
    }
    
    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        unsigned numElements = this->problem_.numElements();
        ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pN =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoN =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobN =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracAinW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracAinN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracWinW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracWinN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *temperature  = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *phaseState   = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityX    = writer.template createField<Scalar, 1>(numElements);
        ScalarField *velocityY    = writer.template createField<Scalar, 1>(numElements);
        ScalarField *velocityZ    = writer.template createField<Scalar, 1>(numElements);

        LocalFunction tmpSol;
        ElementData   elemDat(BoxTraits::ShapeFunctionSetContainer::maxsize);
        VertexData    tmp;

        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();

        for (; elementIt != endit; ++elementIt)
        {
            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            updateElementData_(elemDat, tmpSol, false);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                (*pW)[globalIdx] = elemDat[i].pressure[wPhase];
                (*pN)[globalIdx] = elemDat[i].pressure[nPhase];
                (*pC)[globalIdx] = elemDat[i].pC;
                (*Sw)[globalIdx] = elemDat[i].saturation[wPhase];
                (*Sn)[globalIdx] = elemDat[i].saturation[nPhase];
                (*rhoW)[globalIdx] = elemDat[i].density[wPhase];
                (*rhoN)[globalIdx] = elemDat[i].density[nPhase];
                (*mobW)[globalIdx] = elemDat[i].mobility[wPhase];
                (*mobN)[globalIdx] = elemDat[i].mobility[nPhase];
                (*massfracAinW)[globalIdx] = elemDat[i].massfrac[nComp][wPhase];
                (*massfracAinN)[globalIdx] = elemDat[i].massfrac[nComp][nPhase];
                (*massfracWinW)[globalIdx] = elemDat[i].massfrac[wComp][wPhase];
                (*massfracWinN)[globalIdx] = elemDat[i].massfrac[wComp][nPhase];
                (*temperature)[globalIdx] = asImp_()->temperature((*globalSol)[globalIdx]);
                (*phaseState)[globalIdx] = staticVertexDat_[globalIdx].phaseState;
            };

            // Vector containing the velocity at the element
            GlobalPosition velocity[numPhases];
            GlobalPosition elementVelocity[numPhases];

            // loop over the phases
            for (int phase=0; phase < numPhases; phase++)
            {
                elementVelocity[phase] = 0;

                int elementIdx = this->problem_.elementIdx(*elementIt);
                for (int faceIdx = 0; faceIdx< this->curElementGeom_.numEdges; faceIdx++)
                {
                    velocity[phase] = 0;
                    /* asImp_().calculateDarcyVelocity(velocity[phase],
                                                     faceIdx);
                    */
                    elementVelocity[phase] += velocity[phase];
                }
                elementVelocity[phase] *= 1.0/this->curElementGeom_.numEdges;
                (*velocityX)[elementIdx] = elementVelocity[0][phase];
                if (dim >= 2)
                    (*velocityY)[elementIdx] = elementVelocity[1][phase];
                if (dim == 3)
                    (*velocityZ)[elementIdx] = elementVelocity[2][phase];
            }
        }

        writer.addVertexData(pW, "pW");
        writer.addVertexData(pN, "pN");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(Sw, "SW");
        writer.addVertexData(Sn, "SN");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(rhoN, "rhoN");
        writer.addVertexData(mobW, "mobW");
        writer.addVertexData(mobN, "mobN");
        writer.addVertexData(massfracAinW, "XaW");
        writer.addVertexData(massfracAinN, "XaN");
        writer.addVertexData(massfracWinW, "XwW");
        writer.addVertexData(massfracWinN, "XwN");
        writer.addVertexData(temperature, "T");
        writer.addVertexData(phaseState, "phase state");
        writer.addCellData(velocityX, "Vx");
        if (dim >= 2)
            writer.addCellData(velocityY, "Vy");
        if (dim == 3)
            writer.addCellData(velocityZ, "Vz");
    }

    /*!
     * \brief The storage term of heat
     */
    void heatStorage(SolutionVector &result,
                     int scvIdx,
                     const ElementData &elementData) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Update the temperature gradient at a face of a FV
     *        element.
     */
    void updateTempGrad(GlobalPosition &tempGrad,
                        const GlobalPosition &feGrad,
                        int vertexIdx) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Sets the temperature term of the flux vector to the
     *        heat flux due to advection of the fluids.
     */
    void advectiveHeatFlux(SolutionVector &flux,
                           const PhasesVector &darcyOut,
                           Scalar alpha, // upwind parameter
                           const VertexData *upW, // up/downstream verts
                           const VertexData *dnW,
                           const VertexData *upN,
                           const VertexData *dnN) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     */
    void diffusiveHeatFlux(SolutionVector &flux,
                           int faceIdx,
                           const GlobalPosition &tempGrad) const
    {
        // only relevant for the non-isothermal model!
    }

    /*!
     * \brief Returns the temperature.
     *
     * This is an isothermal model, we let the problem decide which
     * temperature we've got.
     */
    Scalar temperature(const SolutionVector &sol) const
    { return this->problem_.temperature(); }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     */
    void deserializeEntity(std::istream &inStream,
                           const Vertex &vert)
    {
        int vertIdx = this->problem_.vertexIdx(vert);

        // read phase state
        if (!inStream.good()) {
            DUNE_THROW(IOError,
                       "Could not deserialize vertex "
                       << vertIdx);
        }

        inStream >> staticVertexDat_[vertIdx].phaseState;
        staticVertexDat_[vertIdx].oldPhaseState
            = staticVertexDat_[vertIdx].phaseState;
    };

    /*!
     * \brief Write the current phase state of an vertex to a restart
     *        file.
     */
    void serializeEntity(std::ostream &outStream,
                         const Vertex &vert)
    {
        int vertIdx = this->problem_.vertexIdx(vert);

        if (!outStream.good()) {
            DUNE_THROW(IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }
        
        outStream << staticVertexDat_[vertIdx].phaseState
                  << " ";
    };


protected:
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }


    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SpatialFunction &globalSol,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        int phaseState    = staticVertexDat_[globalIdx].phaseState;
        int newPhaseState = phaseState;
        Scalar pW = 0;
        Scalar pN = 0;
        Scalar satW = 0;
        Scalar satN = 0;

        if (formulation == pWsN)
        {
            pW = (*globalSol)[globalIdx][pressureIdx];
            if      (phaseState == bothPhases) satN = (*globalSol)[globalIdx][switchIdx];
            else if (phaseState == wPhaseOnly) satN = 0.0;
            else if (phaseState == nPhaseOnly) satN = 1.0;

                satW = 1 - satN;
            
            LocalPosition local(0);// HACK
            Scalar pC = this->problem_.materialLaw().pC(satW,
                                                        globalPos,
                                                        this->curElement_(), // HACK
                                                        local);
            pN = pW + pC;
        }

        // Evaluate saturation and pressures
        else if (formulation == pNsW)
        {
            pN = (*globalSol)[globalIdx][pressureIdx];
            satW = 0.0;
            if      (phaseState == bothPhases) satW = (*globalSol)[globalIdx][switchIdx];
            else if (phaseState == wPhaseOnly) satW = 1.0;
            else if (phaseState == nPhaseOnly) satW = 0.0;

            satN = 1 - satW;

            LocalPosition local(0);// HACK
            Scalar pC = this->problem_.materialLaw().pC(satW,
                                                        globalPos,
                                                        this->curElement_(), // HACK
                                                        local);
            pW = pN - pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        Scalar temperature = asImp_()->temperature((*globalSol)[globalIdx]);

        // phase state is checked and a switch is performed
        if (phaseState == nPhaseOnly)
        {
            Scalar xWN = (*globalSol)[globalIdx][switchIdx];
            Scalar xWNmax = this->problem_.multicomp().xWN(pN, temperature);

            if (xWN > xWNmax)
            {
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                newPhaseState = bothPhases;
                if (formulation == pNsW) (*globalSol)[globalIdx][switchIdx] = 0.0;
                else if (formulation == pWsN) (*globalSol)[globalIdx][switchIdx] = 1.;
            };
        }
        else if (phaseState == wPhaseOnly)
        {
            Scalar xAW = (*globalSol)[globalIdx][switchIdx];
            Scalar xAWmax = this->problem_.multicomp().xAW(pN, temperature);

            if (xAW > xAWmax)
            {
                // non-wetting phase appears
                std::cout << "Non-wetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                if (formulation == pNsW) (*globalSol)[globalIdx][switchIdx] = 1.0;
                else if (formulation == pWsN) (*globalSol)[globalIdx][switchIdx] = 0.0;
                newPhaseState = bothPhases;
            }
        }
        else if (phaseState == bothPhases) {
            if (satN < 0) {
                // non-wetting phase disappears
                std::cout << "Non-wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                (*globalSol)[globalIdx][switchIdx]
                    = this->problem_.multicomp().xAW(pN, temperature);
                newPhaseState = wPhaseOnly;
            }
            else if (satW < 0) {
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                (*globalSol)[globalIdx][switchIdx]
                    = this->problem_.multicomp().xWN(pN, temperature);
                newPhaseState = nPhaseOnly;
            }
        }

        staticVertexDat_[globalIdx].phaseState = newPhaseState;

        return phaseState != newPhaseState;
    }

    // harmonic mean of the permeability computed directly.  the
    // first parameter is used to store the result.
    static void harmonicMeanK_(Tensor &Ki, const Tensor &Kj)
    {
        for (int kx=0; kx < Tensor::rows; kx++){
            for (int ky=0; ky< Tensor::cols; ky++){
                if (Ki[kx][ky] != Kj[kx][ky]) {
                    Ki[kx][ky] = harmonicMean_(Ki[kx][ky], Kj[kx][ky]);
                }
            }
        }
    }

    // returns the harmonic mean of two scalars
    static Scalar harmonicMean_(Scalar x, Scalar y)
    {
        if (x == 0 || y == 0)
            return 0;
        return (2*x*y)/(x + y);
    };
    
    // parameters given in constructor
    std::vector<StaticVertexData> staticVertexDat_;
    bool                          switchFlag_;
    int                           formulation_;
};


} // end namepace

#endif
