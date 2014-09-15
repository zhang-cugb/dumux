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
#ifndef DUMUX_EVALCFLFLUX_COATS_HH
#define DUMUX_EVALCFLFLUX_COATS_HH

/**
 * @file
 * @brief  Cfl-flux-function to evaluate a Cfl-Condition after Coats 2003
 */

#include <dumux/decoupled/common/impetproperties.hh> 
#include "evalcflflux.hh"

namespace Dumux
{
/*!\ingroup IMPES
 * @brief  Cfl-flux-function to evaluate a Cfl-Condition after Coats 2003
 *
 * tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class EvalCflFluxCoats: public EvalCflFlux<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };
    enum
        {
            wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
            eqIdxPress = Indices::pressureEqIdx,
            eqIdxSat = Indices::satEqIdx,
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
        };

    enum
        {
            pw = Indices::pressureW,
            pn = Indices::pressureNw,
            vt = Indices::velocityTotal,
            sw = Indices::saturationW,
            sn = Indices::saturationNw
        };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

public:
    //! \brief Initializes the cfl-flux-model
    void initialize()
    {
        ElementIterator element = problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(*element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(*element));
        fluidState.setTemperature(problem_.temperature(*element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
    }

    /*! \brief adds a flux to the cfl-criterion evaluation
     *
     * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Element&,int)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux,
                 const Element& element, int phaseIdx = -1)
    {
        addDefaultFlux(flux, phaseIdx);
    }

    /*! \brief adds a flux to the cfl-criterion evaluation
     *
     * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Intersection&,int)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux,
                 const Intersection& intersection, int phaseIdx = -1)
    {
        addDefaultFlux(flux, phaseIdx);
        addCoatsFlux(lambdaW, lambdaNw, viscosityW, viscosityNw, flux, intersection, phaseIdx);
    }

    /*! \brief Returns the Cfl flux-function
     *
     * \copydetails EvalCflFlux::getCflFluxFunction(const Element&)
     */
    Scalar getCflFluxFunction(const Element& element)
    {
    	 Scalar cflFluxDefault = getCflFluxFunctionDefault();

        if (rejectForTimeStepping_)
        	return 0.99 / cflFluxDefault;

        if (std::isnan(cflFluxFunctionCoatsOut_) || std::isinf(cflFluxFunctionCoatsOut_)){cflFluxFunctionCoatsOut_ = 0.0;}
        if (std::isnan(cflFluxFunctionCoatsIn_) || std::isinf(cflFluxFunctionCoatsIn_)){cflFluxFunctionCoatsIn_ = 0.0;}

        Scalar cflFluxFunctionCoats = std::max(cflFluxFunctionCoatsIn_, cflFluxFunctionCoatsOut_);

        if (cflFluxFunctionCoats <= 0)
        {
            return 0.99 / cflFluxDefault;
        }
        else if (cflFluxDefault > cflFluxFunctionCoats)
        {
        	return 0.99 / cflFluxDefault;
        }
        else
        {
            return 0.99 / cflFluxFunctionCoats;
        }
    }

    /*! \brief  Returns the Cfl time-step
     *
     * \copydetails EvalCflFlux::getDt(const Element&)
     */
    Scalar getDt(const Element& element)
    {
        Scalar porosity = std::max(problem_.spatialParams().porosity(element), porosityThreshold_);
        return (getCflFluxFunction(element) * porosity * element.geometry().volume());
    }

    //! \brief  Resets the Timestep-estimator
    void reset()
    {
        cflFluxFunctionCoatsIn_ = 0;
        cflFluxFunctionCoatsOut_ = 0;
        rejectForTimeStepping_ = false;
        fluxWettingOut_ = 0;
        fluxNonwettingOut_ = 0;
        fluxIn_ = 0;
        fluxOut_ = 0;
    }

    /*! \brief Constructs an EvalCflFluxDefault object
     *
     * \param problem A problem type object
     */
    EvalCflFluxCoats(Problem& problem) :
        problem_(problem), epsDerivative_(5e-3), threshold_(1e-8)
    {
        reset();
        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;

        porosityThreshold_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, PorosityThreshold);
    }

private:
    Scalar getCflFluxFunctionDefault()
    {
        if (std::isnan(fluxIn_) || std::isinf(fluxIn_))
        {
            fluxIn_ = 1e-100;
        }

        Scalar cFLFluxIn = fluxIn_;


        Scalar cFLFluxOut = 0;

        if (velocityType_ == vt)
        {
            if (std::isnan(fluxOut_) || std::isinf(fluxOut_))
            {
                fluxOut_ = 1e-100;
            }

            cFLFluxOut = fluxOut_;
        }
        else
        {
            if (std::isnan(fluxWettingOut_) || std::isinf(fluxWettingOut_))
            {
                fluxWettingOut_ = 1e-100;
            }
            if (std::isnan(fluxNonwettingOut_) || std::isinf(fluxNonwettingOut_))
            {
                fluxNonwettingOut_ = 1e-100;
            }

            cFLFluxOut = std::max(fluxWettingOut_, fluxNonwettingOut_);
        }


        //determine timestep
        Scalar cFLFluxFunction = std::max(cFLFluxIn, cFLFluxOut);

        return cFLFluxFunction;
    }

    void addDefaultFlux(Scalar flux,int phaseIdx);

    void addCoatsFlux(Scalar& lambdaW, Scalar& lambdaNw, Scalar& viscosityW, Scalar& viscosityNw, Scalar flux,
                      const Intersection& intersection, int phaseIdx);

    Problem& problem_;//problem data
    Scalar cflFluxFunctionCoatsIn_;
    Scalar cflFluxFunctionCoatsOut_;
    Scalar fluxWettingOut_;
    Scalar fluxNonwettingOut_;
    Scalar fluxOut_;
    Scalar fluxIn_;
    bool rejectForTimeStepping_;
    Scalar density_[numPhases];
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    const Scalar epsDerivative_;
    const Scalar threshold_;
    Scalar porosityThreshold_;
};

template<class TypeTag>
void EvalCflFluxCoats<TypeTag>::addDefaultFlux(Scalar flux, int phaseIdx)
{
    switch (phaseIdx)
    {
    case wPhaseIdx:
        {
            //for time step criterion
            if (flux >= 0)
            {
                fluxWettingOut_ += flux;
            }
            if (flux < 0)
            {
                fluxIn_ -= flux;
            }

            break;
        }

        //for time step criterion if the non-wetting phase velocity is used
    case nPhaseIdx:
        {
            if (flux >= 0)
            {
                fluxNonwettingOut_ += flux;
            }
            if (flux < 0)
            {
                fluxIn_ -= flux;
            }

            break;
        }
    default:
        {
            if (flux >= 0)
            {
                fluxOut_ += flux;
            }
            if (flux < 0)
            {
                fluxIn_ -= flux;
            }

            break;
        }
    }
}

/*! \brief adds a flux to the cfl-criterion evaluation
 *
 * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Intersection&,int)
 */
template<class TypeTag>
void EvalCflFluxCoats<TypeTag>::addCoatsFlux(Scalar& lambdaW, Scalar& lambdaNw,
                                             Scalar& viscosityW, Scalar& viscosityNw, Scalar flux,
                                             const Intersection& intersection, int phaseIdx)
{
    if (rejectForTimeStepping_)
        return;
    else if (phaseIdx != wPhaseIdx)
        return;

    Scalar lambdaT = (lambdaW + lambdaNw);

    if (lambdaT <= threshold_)
        return;

    ElementPointer element = intersection.inside();

    //coordinates of cell center
    const GlobalPosition& globalPos = element->geometry().center();

    // cell index
    int globalIdxI = problem_.variables().index(*element);

    CellData& cellDataI = problem_.variables().cellData(globalIdxI);

    if (cellDataI.potential(wPhaseIdx) < 0.0 || cellDataI.potential(nPhaseIdx) < 0.0)
    {
        rejectForTimeStepping_ = true;
        cflFluxFunctionCoatsIn_ = 0;
        cflFluxFunctionCoatsOut_ = 0;
        return;
    }

    int indexInInside = intersection.indexInInside();

    Scalar satI = cellDataI.saturation(wPhaseIdx);
    Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
    Scalar lambdaNwI = cellDataI.mobility(nPhaseIdx);

    Scalar dpc_dsI = MaterialLaw::dpc_dsw(problem_.spatialParams().materialLawParams(*element), satI);

    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();

    if (intersection.neighbor())
    {
        ElementPointer neighbor = intersection.outside();

        int globalIdxJ = problem_.variables().index(*neighbor);

        CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

        if (cellDataJ.potential(wPhaseIdx) < 0.0 || cellDataJ.potential(nPhaseIdx) < 0.0 )
        {
            rejectForTimeStepping_ = true;
            cflFluxFunctionCoatsIn_ = 0;
            cflFluxFunctionCoatsOut_ = 0;
            return;
        }

        if (element->level() != neighbor->level())
        {
            rejectForTimeStepping_ = true;
            cflFluxFunctionCoatsIn_ = 0;
            cflFluxFunctionCoatsOut_ = 0;
            return;
        }

        bool takeNeighbor = (element->level() < neighbor->level());
        //get phase potentials
        bool upwindWI =
            (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(wPhaseIdx, intersection.indexInOutside()) :
            cellDataI.fluxData().isUpwindCell(wPhaseIdx, indexInInside);
        bool upwindNwI =
            (takeNeighbor) ? !cellDataJ.fluxData().isUpwindCell(nPhaseIdx, intersection.indexInOutside()) :
            cellDataI.fluxData().isUpwindCell(nPhaseIdx, indexInInside);

            const GlobalPosition& globalPosNeighbor = neighbor->geometry().center();

            // distance vector between barycenters
            GlobalPosition distVec = globalPosNeighbor - globalPos;

            // compute distance between cell centers
            Scalar dist = std::abs(distVec * unitOuterNormal);

            Scalar satJ = cellDataJ.saturation(wPhaseIdx);
            Scalar lambdaWJ = cellDataI.mobility(wPhaseIdx);
            Scalar lambdaNwJ = cellDataI.mobility(nPhaseIdx);

            Scalar dpc_dsJ = MaterialLaw::dpc_dsw(problem_.spatialParams().materialLawParams(*neighbor), satJ);

            // compute vectorized permeabilities
            DimVector permeability(0);
            DimMatrix perm(0);
            problem_.spatialParams().meanK(perm, problem_.spatialParams().intrinsicPermeability(*element));
            perm.mv(unitOuterNormal, permeability);

            Scalar perm1 = permeability * unitOuterNormal;

            permeability = 0;
            problem_.spatialParams().meanK(perm, problem_.spatialParams().intrinsicPermeability(*neighbor));
            perm.mv(unitOuterNormal, permeability);

            Scalar perm2 = permeability * unitOuterNormal;

            Scalar meanPermeability = 2*perm1*perm2/(perm1 + perm2);

            Scalar transmissibility =  meanPermeability * intersection.geometry().volume() / dist;

            Scalar satUpw = 0;
            if (upwindWI)
            {
                satUpw = std::max(satI, 0.0);
            }
            else
            {
                satUpw = std::max(satJ, 0.0);
            }

            Scalar ds = epsDerivative_;

            Scalar satPlus = satUpw + epsDerivative_;
            Scalar satMinus = satUpw;
            if (satMinus - epsDerivative_ > 0.0)
            {
                satMinus -= epsDerivative_;
                ds += epsDerivative_;
            }

            Scalar dLambdaWDs = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*neighbor), std::abs(satPlus)) / viscosityW;
            dLambdaWDs -= MaterialLaw::krw(problem_.spatialParams().materialLawParams(*neighbor), std::abs(satMinus)) / viscosityW;
            dLambdaWDs /= (ds);

            if (upwindNwI)
            {
                satUpw = std::max(1 - satI, 0.0);
            }
            else
            {
                satUpw = std::max(1 - satJ, 0.0);
            }

            ds = epsDerivative_;

            satPlus = satUpw + epsDerivative_;
            satMinus = satUpw;
            if (satMinus - epsDerivative_ > 0.0)
            {
                satMinus -= epsDerivative_;
                ds += epsDerivative_;
            }

            Scalar dLambdaNwDs = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*neighbor), satPlus) / viscosityNw;
            dLambdaNwDs -= MaterialLaw::krn(problem_.spatialParams().materialLawParams(*neighbor), satMinus) / viscosityNw;
            dLambdaNwDs /= (ds);

            Scalar lambdaWCap = 0.5 * (lambdaWI + lambdaWJ);
            Scalar lambdaNwCap = 0.5 * (lambdaNwI + lambdaNwJ);

            Scalar potentialDiff = cellDataI.potential(wPhaseIdx) - cellDataJ.potential(wPhaseIdx);
            Scalar cflFlux = transmissibility * lambdaNw * dLambdaWDs * std::abs(potentialDiff) / lambdaT;

            potentialDiff = cellDataI.potential(nPhaseIdx) - cellDataJ.potential(nPhaseIdx);
            cflFlux -= transmissibility * lambdaW * dLambdaNwDs * std::abs(potentialDiff) / lambdaT;

            cflFlux -= transmissibility * lambdaWCap * lambdaNwCap * (dpc_dsI + dpc_dsJ) / lambdaT;

            if ((upwindWI && lambdaW > threshold_)|| (upwindNwI && lambdaW < threshold_))
            {
                cflFluxFunctionCoatsOut_ += cflFlux;
            }
            else
            {
            	cflFluxFunctionCoatsIn_ += cflFlux;
            }
    }
    else
    {
            // center of face in global coordinates
            GlobalPosition globalPosFace = intersection.geometry().center();

            //get boundary type
            BoundaryTypes bcType;
            problem_.boundaryTypes(bcType, intersection);

            // distance vector between barycenters
            Dune::FieldVector < Scalar, dimWorld > distVec = globalPosFace - globalPos;

            // compute distance between cell centers
            Scalar dist = distVec.two_norm();

            //permeability vector at boundary
            // compute vectorized permeabilities

            Dune::FieldVector<Scalar, dim> permeability(0);
            DimMatrix perm(0);
            problem_.spatialParams().meanK(perm, problem_.spatialParams().intrinsicPermeability(*element));
            perm.mv(unitOuterNormal, permeability);          
            
        	Scalar faceArea = intersection.geometry().volume();

        	Scalar transmissibility = (unitOuterNormal * permeability) * faceArea / dist;

            Scalar satWBound =  cellDataI.saturation(wPhaseIdx);
            if (bcType.isDirichlet(eqIdxSat))
            {
                PrimaryVariables bcValues;
                problem_.dirichlet(bcValues, intersection);
                switch (saturationType_)
                {
                case sw:
                    {
                        satWBound = bcValues[eqIdxSat];
                        break;
                    }
                case sn:
                    {
                        satWBound = 1 - bcValues[eqIdxSat];
                        break;
                    }
                default:
                    {
                        DUNE_THROW(Dune::RangeError, "saturation type not implemented");
                        break;
                    }
                }

            }
            
            Scalar potWBound =  cellDataI.potential(wPhaseIdx);
        	Scalar potNwBound =  cellDataI.potential(nPhaseIdx);
        	Scalar gdeltaZ = (problem_.bBoxMax()-globalPosFace) * problem_.gravity();
        	if (bcType.isDirichlet(eqIdxPress))
        	{
            	PrimaryVariables bcValues;
            	problem_.dirichlet(bcValues, intersection);
            	switch (pressureType_)
            	{
            	case pw:
                	{
                    	potWBound = bcValues[eqIdxPress] + density_[wPhaseIdx] * gdeltaZ;
                       potNwBound = bcValues[eqIdxPress] + MaterialLaw::pc(problem_.spatialParams().materialLawParams(*element), satWBound)
                        								 + density_[nPhaseIdx] * gdeltaZ;
                    	break;
                	}
            	case pn:
                	{
                           potWBound = bcValues[eqIdxPress] - MaterialLaw::pc(problem_.spatialParams().materialLawParams(*element),satWBound)
                           									+ density_[wPhaseIdx] * gdeltaZ;
                    	potNwBound = bcValues[eqIdxPress] + density_[nPhaseIdx] * gdeltaZ;
                   		break;
                	}
            	default:
                	{
                	    DUNE_THROW(Dune::RangeError, "pressure type not implemented");
                    	break;
                	}
            	}
        	}
        	else if (bcType.isNeumann(eqIdxPress) && bcType.isDirichlet(eqIdxSat))
        	{
            	PrimaryVariables bcValues;
            	problem_.neumann(bcValues, intersection);

	            bcValues[wPhaseIdx] /= density_[wPhaseIdx];
    	        bcValues[nPhaseIdx] /= density_[nPhaseIdx];
	
    	        bcValues[wPhaseIdx] *= faceArea;
        	    bcValues[nPhaseIdx] *= faceArea;
	
    	        bool hasPotWBound = false;
        	    if (lambdaW != 0 && bcValues[wPhaseIdx] != 0)
            	{
        	        potWBound -= bcValues[wPhaseIdx] / (transmissibility * lambdaW);
        	        hasPotWBound = true;
        	    }
      	     	bool hasPotNwBound = false;
        	    if (lambdaNw != 0 && bcValues[nPhaseIdx] != 0)
            	{
      	          	potNwBound -= bcValues[nPhaseIdx] / (transmissibility * lambdaNw);
        	        hasPotNwBound = true;
            	}

  	          	if (hasPotWBound && !hasPotNwBound)
    	      	{
                       potNwBound = potWBound + MaterialLaw::pc(problem_.spatialParams().materialLawParams(*element),satWBound)
                       						  + (density_[nPhaseIdx] - density_[wPhaseIdx]) * gdeltaZ;
   	          	}
    	        else if (!hasPotWBound && hasPotNwBound)
        	    {
                   potWBound = potNwBound - MaterialLaw::pc(problem_.spatialParams().materialLawParams(*element),satWBound)
                   						  + (density_[nPhaseIdx] - density_[wPhaseIdx]) * gdeltaZ;
            	}
        	}
            else if (bcType.isNeumann(eqIdxPress))
            {
                PrimaryVariables bcValues;
                problem_.neumann(bcValues, intersection);

                bcValues[wPhaseIdx] /= density_[wPhaseIdx];
                bcValues[nPhaseIdx] /= density_[nPhaseIdx];

                bcValues[wPhaseIdx] *= faceArea;
                bcValues[nPhaseIdx] *= faceArea;

                if (bcValues[wPhaseIdx] > 0)
                {
                    cflFluxFunctionCoatsOut_ += std::abs(bcValues[wPhaseIdx]);
                }
                else
                {
                    cflFluxFunctionCoatsIn_ += std::abs(bcValues[wPhaseIdx]);
                }
                if (bcValues[nPhaseIdx] > 0)
                {
                    cflFluxFunctionCoatsOut_ += std::abs(bcValues[nPhaseIdx]);
                }
                else
                {
                    cflFluxFunctionCoatsIn_ += std::abs(bcValues[nPhaseIdx]);
                }

                return;
            }
            else
            {
                rejectForTimeStepping_ = true;
                cflFluxFunctionCoatsIn_ = 0;
                cflFluxFunctionCoatsOut_ = 0;
                return;
            }

            Scalar dpc_dsBound = MaterialLaw::dpc_dsw(problem_.spatialParams().materialLawParams(*element), satWBound);

            Scalar lambdaWBound = 0;
            Scalar lambdaNwBound = 0;

            Scalar temperature = problem_.temperature(*element);
            Scalar referencePressure = problem_.referencePressure(*element);
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, referencePressure);
            fluidState.setPressure(nPhaseIdx, referencePressure);
            fluidState.setTemperature(temperature);

            Scalar viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx);
            Scalar viscosityNwBound =
                FluidSystem::viscosity(fluidState, nPhaseIdx);
            lambdaWBound = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satWBound) / viscosityWBound;
            lambdaNwBound = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satWBound) / viscosityNwBound;

            Scalar satUpw = 0;
            if (cellDataI.fluxData().isUpwindCell(wPhaseIdx, indexInInside))
            {
                satUpw = std::max(satI, 0.0);
            }
            else
            {
                satUpw = std::max(satWBound, 0.0);
            }

            Scalar ds = epsDerivative_;

            Scalar satPlus = satUpw + epsDerivative_;
            Scalar satMinus = satUpw;
            if (satMinus - epsDerivative_ > 0.0)
            {
                satMinus -= epsDerivative_;
                ds += epsDerivative_;
            }

            Scalar dLambdaWDs = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satPlus) / viscosityW;
            dLambdaWDs -= MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satMinus) / viscosityW;
            dLambdaWDs /= (ds);

            if (cellDataI.fluxData().isUpwindCell(nPhaseIdx, indexInInside))
            {
                satUpw = std::max(1 - satI, 0.0);
            }
            else
            {
                satUpw = std::max(1 - satWBound, 0.0);
            }

            ds = epsDerivative_;

            satPlus = satUpw + epsDerivative_;
            satMinus = satUpw;
            if (satMinus - epsDerivative_ > 0.0)
            {
                satMinus -= epsDerivative_;
                ds += epsDerivative_;
            }

            Scalar dLambdaNwDs = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satPlus) / viscosityNw;
            dLambdaNwDs -= MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satMinus) / viscosityNw;
            dLambdaNwDs /= (ds);

            Scalar lambdaWCap = 0.5 * (lambdaWI + lambdaWBound);
            Scalar lambdaNwCap = 0.5 * (lambdaNwI + lambdaNwBound);

            Scalar potDiff = cellDataI.potential(wPhaseIdx) - potWBound;
            Scalar cflFlux = transmissibility * lambdaNw * dLambdaWDs * std::abs(potDiff) / lambdaT;

            cflFlux -= transmissibility * lambdaWCap * lambdaNwCap * (dpc_dsI + dpc_dsBound) / lambdaT;

            potDiff = cellDataI.potential(nPhaseIdx) - potNwBound;
            cflFlux -= transmissibility * lambdaW * dLambdaNwDs * std::abs(potDiff) / lambdaT;

            if ((cellDataI.fluxData().isUpwindCell(wPhaseIdx, indexInInside) && lambdaW > threshold_) ||
                (cellDataI.fluxData().isUpwindCell(nPhaseIdx, indexInInside) && lambdaW < threshold_))
            {
                cflFluxFunctionCoatsOut_ += cflFlux;
            }
            else
            {
            	cflFluxFunctionCoatsIn_ += cflFlux;
            }
    }
}

}

#endif
