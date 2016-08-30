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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase one-component fully implicit model.
 *
 *Important note: The 2p1c model requires the use of the non-isothermal extension found in dumux/implicit/nonisothermal
 */
#ifndef DUMUX_2P1C_LOCAL_RESIDUAL_HH
#define DUMUX_2P1C_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPOneCNIModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase one-component fully implicit model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 2P1C flow.
 */
template<class TypeTag>
class TwoPOneCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;


    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        conti0EqIdx = Indices::conti0EqIdx,//!< Index of the mass conservation equation for the water component
        energyEqIdx = Indices::energyEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    static const bool useBlockingOfSpuriousFlow = GET_PROP_VALUE(TypeTag, UseBlockingOfSpuriousFlow);

public:

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        storage = 0;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                storage[conti0EqIdx] +=
                    volVars.porosity()
                    * volVars.saturation(phaseIdx) * volVars.density(phaseIdx);
            }
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param faceIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, const int faceIdx, const bool onBoundary=false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
                           this->element_(),
                           this->fvGeometry_(),
                           faceIdx,
                           this->curVolVars_(),
                           onBoundary);

        flux = 0;
        asImp_()->computeMyAdvectiveFlux(flux, fluxVars); //Method is called "computeMyAdvectiveFlux"
                                                          // to prevent overwrite by nilocalresidual.hh
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */

    void computeMyAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        Scalar massUpwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

            Scalar factor = 1.0;

            // factor for the blocking of spurious cold water fluxes into the steam zone (Gudbjerg, 2004)
            if(useBlockingOfSpuriousFlow)
            { factor = factor_(up, dn, phaseIdx); }


                // add advective flux of current component in current
                // phase
                // if alpha > 0 und alpha < 1 then both upstream and downstream
                // nodes need their contribution
                // if alpha == 1 (which is mostly the case) then, the downstream
                // node is not evaluated
                int eqIdx = conti0EqIdx;
                flux[eqIdx] +=  fluxVars.volumeFlux(phaseIdx)
                        * (massUpwindWeight * up.fluidState().density(phaseIdx)
                        + (1.0 - massUpwindWeight) * dn.fluidState().density(phaseIdx) )
                        * factor ;
        }

        // advective heat flux in all phases
        flux[energyEqIdx] = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // vertex data of the upstream and the downstream vertices
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

            Scalar factor = 1.0;

            if(useBlockingOfSpuriousFlow)
            { factor = factor_(up, dn, phaseIdx); }

            flux[energyEqIdx] +=
                fluxVars.volumeFlux(phaseIdx) * (
                                              massUpwindWeight * // upstream vertex
                                              up.density(phaseIdx) *
                                              up.enthalpy(phaseIdx)
                                              +
                                              (1-massUpwindWeight) * // downstream vertex
                                              dn.density(phaseIdx) *
                                              dn.enthalpy(phaseIdx) )
                                              * factor ;
        }
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {}

    const bool spuriousFlowDetected() const
    {
        return spuriousFlowDetected_;
    }

    void resetSpuriousFlowDetected()
    {
        spuriousFlowDetected_ = false;
    }

protected:

    void evalPhaseStorage_(const int phaseIdx)
    {
        // evaluate the storage terms of a single phase
        for (int i=0; i < this->fvGeometry_().numScv; i++) {
            PrimaryVariables &storage = this->storageTerm_[i];
            const ElementVolumeVariables &elemVolVars = this->curVolVars_();
             const VolumeVariables &volVars = elemVolVars[i];

            // compute storage term of all components within all phases
           storage = 0;
           storage[conti0EqIdx] += volVars.density(phaseIdx)
                        * volVars.saturation(phaseIdx);

           storage *= volVars.porosity();
           storage *= this->fvGeometry_().subContVol[i].volume;
          }
    }

    /*!
     * \brief Calculate the blocking factor which prevents spurious cold water fluxes into the steam zone (Gudbjerg, 2005)
     *
     * \param up The upstream volume variables
     * \param dn The downstream volume variables
     * \param phaseIdx The index of the fluid phase
     */
    Scalar factor_(const VolumeVariables &up, const VolumeVariables &dn,const  int phaseIdx) const {

        Scalar factor = 1.0;

        Scalar tDn = dn.temperature(); //temperature of the downstream SCV (where the cold water is potentially intruding into a steam zone)
        Scalar tUp = up.temperature(); //temperature of the upstream SCV

        Scalar sgDn = dn.saturation(gPhaseIdx); //gas phase saturation of the downstream SCV
        Scalar sgUp = up.saturation(gPhaseIdx); //gas phase saturation of the upstream SCV

        bool upIsNotSteam = false;
        bool downIsSteam = false;
        bool spuriousFlow = false;

        if(sgUp <= 1e-5)
        {upIsNotSteam = true;}

        if(sgDn > 1e-5)
        {downIsSteam = true;}

        if(upIsNotSteam && downIsSteam  && tDn > tUp && phaseIdx == wPhaseIdx)
        {
          spuriousFlow = true;
          spuriousFlowDetected_ = true;
        }

        if(spuriousFlow)
        {
          Scalar deltaT = tDn - tUp;

          if((deltaT) > 15 )
          {factor = 0.0 ;}

          else
          { factor = 1-(deltaT/15); }

        }
        return factor;
    }

    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }

    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }

private:
    mutable bool spuriousFlowDetected_;// true if spurious cold-water flow into a steam zone is happening
};

} // end namespace

#endif
