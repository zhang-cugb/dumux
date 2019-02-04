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
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
#ifndef DUMUX_INJECTION_SPATIAL_PARAMETERS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMETERS_HH

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotmateriallaw.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/porosityprecipitation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCMinTests
 * \brief Spatial parameters for the dissolution problem
 * where water is injected in a for flushing precipitated salt clogging a gas reservoir.
 */
template<class FVGridGeometry, class Scalar>
class DissolutionSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         DissolutionSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       DissolutionSpatialParams<FVGridGeometry, Scalar>>;

    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    //! export the material law type used
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

    DissolutionSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        solubilityLimit_       = getParam<Scalar>("SpatialParams.SolubilityLimit", 0.26);
        referencePorosity_     = getParam<Scalar>("SpatialParams.referencePorosity", 0.11);
        referencePermeabilityLeft_ = getParam<Scalar>("SpatialParams.referencePermeabilityLeft", 2.23e-14);
        referencePermeabilityRight_ = getParam<Scalar>("SpatialParams.referencePermeabilityRight", 2.23e-14);
        irreducibleLiqSat_     = getParam<Scalar>("SpatialParams.IrreducibleLiqSat", 0.2);
        irreducibleGasSat_     = getParam<Scalar>("SpatialParams.IrreducibleGasSat", 1e-3);
        vgAlphaLeft_               = getParam<Scalar>("SpatialParams.VGAlphaLeft", 1.5);
        vgnLeft_                   = getParam<Scalar>("SpatialParams.VGnLeft", 4.0);
        vgAlphaRight_               = getParam<Scalar>("SpatialParams.VGAlphaRight", 1.5);
        vgnRight_                   = getParam<Scalar>("SpatialParams.VGnRight", 4.0);

        plotFluidMatrixInteractions_ = getParam<bool>("Output.PlotFluidMatrixInteractions");

        //Van Genuchen parameters
        // residual saturations
        materialParams_.setSwr(irreducibleLiqSat_);
        materialParams_.setSnr(irreducibleGasSat_);

        //Van Genuchen parameters
        materialParams_.setVgAlpha(vgAlphaLeft_ );
        materialParams_.setVgn(vgnLeft_);

        // residual saturations
        materialParamsCoarse_.setSwr(irreducibleLiqSat_);
        materialParamsCoarse_.setSnr(irreducibleGasSat_);

        //Van Genuchen parameters
        materialParamsCoarse_.setVgAlpha(vgAlphaRight_);
        materialParamsCoarse_.setVgn(vgnLeft_);
                //Van Genuchen parameters
//         materialParamsCoarse_.setVgAlpha(vgAlpha_ );
//         materialParamsCoarse_.setVgn(vgn_);
    }

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        PlotMaterialLaw<Scalar, MaterialLaw> plotMaterialLaw;
        GnuplotInterface<Scalar> gnuplot(plotFluidMatrixInteractions_);
        gnuplot.setOpenPlotWindow(plotFluidMatrixInteractions_);
        plotMaterialLaw.addpcswcurve(gnuplot, materialParams_, irreducibleLiqSat_, 1-irreducibleGasSat_);
        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set yrange [0:15000]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("pc-Sw");

        gnuplot.resetAll();
        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set yrange [0:1]");
        plotMaterialLaw.addkrcurves(gnuplot, materialParams_, irreducibleLiqSat_, 1/*-irreducibleGasSat_*/, "fine");
        gnuplot.plot("kr");
    }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minimalPorosity(const Element& element, const SubControlVolume &scv) const
    { return 1e-5; }

    /*!
     *  \brief Define the volume fraction of the inert component
     *
     *  \param globalPos The global position in the domain
     *  \param compIdx The index of the inert solid component
     */
    template<class SolidSystem>
    Scalar inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { return 1.0-referencePorosity_; }

    /*!
     *  \brief Define the reference porosity \f$[-]\f$ distribution.
     *  This is the porosity of the porous medium without any of the
     *  considered solid phases.
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    {
         Scalar pos0 = scv.center()[0];
         if (pos0<0.035+eps_){
            return referencePorosity_;
         }
         else
             return referencePorosity_;
        }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *
     *  Solution dependent permeability function
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center());

        Scalar sumPrecipitates = 0.0;
        sumPrecipitates += priVars[3 /*numComp*/];

         using std::max;
         const auto poro =  max(/*minPoro*/1e-5, referencePorosity_ - sumPrecipitates);

         Scalar pos0 = scv.center()[0];
         if (pos0<0.035+eps_){
         return permLaw_.evaluatePermeability(referencePermeabilityLeft_, referencePorosity_, poro);
         }
         else
         {
             return permLaw_.evaluatePermeability(referencePermeabilityRight_, referencePorosity_, poro);
         }
    }

//     Scalar solidity(const SubControlVolume &scv) const
//     { return 1.0 - porosityAtPos(scv.center()); }

    Scalar solubilityLimit() const
    { return solubilityLimit_; }

    Scalar theta(const SubControlVolume &scv) const
    { return 10.0; }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0]<0.035+eps_){
            return materialParams_;
        }
        else
        {
            return materialParamsCoarse_;
        }

    }

    // define which phase is to be considered as the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:

    MaterialLawParams materialParams_;
    MaterialLawParams materialParamsCoarse_;

    PermeabilityKozenyCarman<PermeabilityType> permLaw_;

    Scalar solubilityLimit_;
    Scalar referencePorosity_;
    PermeabilityType referencePermeabilityLeft_ = 0.0;
    PermeabilityType referencePermeabilityRight_ = 0.0;
    Scalar irreducibleLiqSat_;
    Scalar irreducibleGasSat_;
    Scalar vgAlphaLeft_;
    Scalar vgnLeft_ ;
    Scalar vgAlphaRight_;
    Scalar vgnRight_ ;
    static constexpr Scalar eps_ = 1e-6;

    bool plotFluidMatrixInteractions_;
};

} // end namespace Dumux

#endif
