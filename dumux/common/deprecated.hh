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
 * \brief Helpers for deprecation
 */

#ifndef DUMUX_COMMON_DEPRECATED_HH
#define DUMUX_COMMON_DEPRECATED_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {

#ifndef DOXYGEN // hide from doxygen
// Helper classes/functions for deprecation
// Each implementation has to state after which release
// it will be removed. Implementations in the Deprecated
// namespace will be removed without
// deprecation after their usage in the code exprired,
// so most likely you don't want to use this in your code
namespace Deprecated {

////////////////////////////////////////////
// Remove the following after Release 3.2 //
////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// Deprecation warnings for effective diffusion coefficients //
///////////////////////////////////////////////////////////////
constexpr auto hasEffDiffCoeffImpl = Dumux::isValid([](auto&& v) -> decltype(v.effectiveDiffusionCoefficient(0,0,0)){return 0;});
template<class VolumeVariables> constexpr bool hasEffDiffCoeff = decltype(hasEffDiffCoeffImpl(std::declval<VolumeVariables>())){};

template<class EDL, class VV,
         typename std::enable_if_t<!hasEffDiffCoeff<VV>, int> = 0>
[[deprecated("The volume variables class used does not have an effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx) function. "
             "This will become mandatory after the 3.2 release! See e.g. the files /dumux/porousmediumflow/2p2c/volumevariables.hh "
             "and /dumux/porousmediumflow/2p2c/properties.hh to see how this can be realized.")]]
decltype(EDL::effectiveDiffusionCoefficient(std::declval<VV>(), int(), int(), int()))
effectiveDiffusionCoefficient(const VV& volVars, int phaseIdx, int compIIdx, int compJIdx)
{
    return EDL::effectiveDiffusionCoefficient(volVars, phaseIdx, compIIdx, compJIdx);
}

template<class EDL, class VV,
    typename std::enable_if_t<hasEffDiffCoeff<VV>, int> = 0>
auto effectiveDiffusionCoefficient(const VV& volVars, int phaseIdx, int compIIdx, int compJIdx)
{ return volVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx); }

    ////////////////////
    // Maxwell Stefan //
    ////////////////////

template<class EDL, class FS, class VV, class P, class E, class SCV,
    typename std::enable_if_t<hasEffDiffCoeff<VV>, int> = 0>
auto effectiveMSDiffusionCoefficient(const VV& volVars,
                                     const int phaseIdx,
                                     const int compIIdx,
                                     const int compJIdx,
                                     const P& problem,
                                     const E& element,
                                     const SCV& scv)
{ return volVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx); }

template<class EDL, class FS, class VV, class P, class E, class SCV,
    typename std::enable_if_t<!hasEffDiffCoeff<VV>, int> = 0>
[[deprecated("The volume variables class used does not have an effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx) function. "
             "This will become mandatory after the 3.2 release! See e.g. the files /dumux/porousmediumflow/2p2c/volumevariables.hh "
             "and /dumux/porousmediumflow/2p2c/properties.hh to see how this can be realized.")]]
auto effectiveMSDiffusionCoefficient(const VV& volVars,
                                     const int phaseIdx,
                                     const int compIIdx,
                                     const int compJIdx,
                                     const P& problem,
                                     const E& element,
                                     const SCV& scv)
{
    auto tinInside = 0.0;
    if constexpr (FS::isTracerFluidSystem())
    {
        tinInside = FS::binaryDiffusionCoefficient(compIIdx,
                                                   compJIdx,
                                                   problem,
                                                   element,
                                                   scv);
    }
    else
    {
        auto fluidState = volVars.fluidState();
        typename FS::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        tinInside = FS::binaryDiffusionCoefficient(fluidState,
                                                   paramCache,
                                                   phaseIdx,
                                                   compIIdx,
                                                   compJIdx);
    }

    return EDL::effectiveDiffusivity(volVars.porosity(), volVars.saturation(phaseIdx), tinInside);
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

} // end namespace Deprecated
#endif

} // end namespace Dumux

#endif
