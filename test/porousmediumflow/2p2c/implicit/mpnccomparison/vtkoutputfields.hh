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
 * \ingroup TwoPTwoCModel
 * \brief Adds vtk output fields specific to the twop model
 */
#ifndef DUMUX_TWOPTWOC_MPNC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_TWOPTWOC_MPNC_VTK_OUTPUT_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Adds vtk output fields specific to the two-phase two-component model
 */
class TwoPTwoCMPNCVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized vtk output fields
        vtk.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
        vtk.addVolumeVariable([](const auto& v){ return v.saturation(FluidSystem::phase1Idx); }, "Sn");
        vtk.addVolumeVariable([](const auto& v){ return v.saturation(FluidSystem::phase0Idx); }, "Sw");
        vtk.addVolumeVariable([](const auto& v){ return v.pressure(FluidSystem::phase1Idx); }, "pn");
        vtk.addVolumeVariable([](const auto& v){ return v.pressure(FluidSystem::phase0Idx); }, "pw");

        vtk.addVolumeVariable([](const auto& v){ return v.density(FluidSystem::phase0Idx); }, "rhoW");
        vtk.addVolumeVariable([](const auto& v){ return v.density(FluidSystem::phase1Idx); }, "rhoN");
        vtk.addVolumeVariable([](const auto& v){ return v.mobility(FluidSystem::phase0Idx); }, "mobW");
        vtk.addVolumeVariable([](const auto& v){ return v.mobility(FluidSystem::phase1Idx); }, "mobN");

        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                vtk.addVolumeVariable([i,j](const auto& v){ return v.massFraction(i,j); },"X_"+ FluidSystem::phaseName(i) + "^" + FluidSystem::componentName(j));

        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                vtk.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },"x_"+ FluidSystem::phaseName(i) + "^" + FluidSystem::componentName(j));
    }
};

} // end namespace Dumux

#endif
