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
 * \brief Class to store time level information during time integration.
 * \todo TODO: Properly implement this thing and give it a good name
 *             Also, this should probably go to the common folder
 */
#ifndef DUMUX_TIME_LEVEL_HH
#define DUMUX_TIME_LEVEL_HH

namespace Dumux {

template<class Scalar>
class TimeLevel
{
public:

    TimeLevel(Scalar curTime)
    : curTime_(curTime)
    , prevTime_(curTime)
    , timeStepFraction_(1.0)
    {}

    TimeLevel(Scalar curTime,
              Scalar prevTime,
              Scalar dtFraction)
    : curTime_(curTime)
    , prevTime_(prevTime)
    , timeStepFraction_(dtFraction)
    {}

    Scalar current() const
    { return curTime_; }

    Scalar previous() const
    { return prevTime_; }

    Scalar timeStepFraction() const
    { return timeStepFraction_; }

private:
    Scalar curTime_;
    Scalar prevTime_;
    Scalar timeStepFraction_;
};

} // end namespace Dumux

#endif
