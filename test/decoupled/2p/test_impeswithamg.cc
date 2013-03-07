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
 * \ingroup IMPETtests
 * \brief test for the sequential 2p model
 */
#include "config.h"

#if HAVE_DUNE_PDELAB

// Check if DUNE-PDELab has been patched for our needs. 
// TODO: this should be replaced by a proper configure test.
#include <dune/pdelab/backend/istlvectorbackend.hh>
#ifdef DUNE_PDELAB_IS_PATCHED_FOR_DUMUX

#include "test_impesproblem.hh"
#include <dumux/common/start.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe list of mandatory arguments for this program is:\n"
                                       "\t-TimeManager.TEnd      End of the simulation [s] \n"
                                       "\t-TimeManager.DtInitial Initial timestep size [s] \n"
                                       "\t-Grid.NumberOfCellsX   Resolution in x-direction [-]\n"
                                       "\t-Grid.NumberOfCellsY   Resolution in y-direction [-]\n"
                                       "\t-Grid.UpperRightX      Length of the domain [m]\n"
                                       "\t-Grid.UpperRightY      Height of the domain [m]\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    typedef TTAG(IMPESTestProblemWithAMG) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, usage);
}
#else // DUNE_PDELAB_IS_PATCHED_FOR_DUMUX

#warning You need to have  a patched dune-pdelab to run this test, see ../../../patches/README for details.

#include <iostream>

int main()
{
    std::cerr << "You need to have a patched dune-pdelab to run this test, "
                 "see ../../../patches/README for details." << std::endl;;
    return 77;
}

#endif // DUNE_PDELAB_IS_PATCHED_FOR_DUMUX

#else // HAVE_DUNE_PDELAB

#warning You need to have dune-pdelab installed and patched to run this test.

#include <iostream>

int main()
{
    std::cerr << "You need to have dune-pdelab installed and patched to run this test.\n";
    return 77;
}
#endif // HAVE_DUNE_PDELAB
