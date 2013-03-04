/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 20010 by Markus Wolff                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Test for the adaptive 2p2c model
 */
#include "config.h"

#include "test_adaptive2p2cproblem.hh"
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
                    errorMessageOut += "\n\nThe list of mandatory options for this program is:\n"
                                        "\t-TimeManager.TEnd              End of the simulation [s] \n"
                                        "\t-TimeManager.DtInitial         Initial timestep size [s] \n"
                                        "\t-Grid.File                     Name of the file containing the grid \n"
                                        "\t                               definition in DGF format\n"
                                        "\t-SpatialParams.LensLowerLeftX  x-coordinate of the lower left corner of the lens [m] \n"
                                        "\t-SpatialParams.LensLowerLeftY  y-coordinate of the lower left corner of the lens [m] \n"
                                        "\t-SpatialParams.LensUpperRightX x-coordinate of the upper right corner of the lens [m] \n"
                                        "\t-SpatialParams.LensUpperRightY y-coordinate of the upper right corner of the lens [m] \n"
                                        "\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}
// The main function using the standard start procedure
int main(int argc, char** argv)
{
#if HAVE_ALUGRID
    typedef TTAG(Adaptive2p2c) TypeTag;
    return Dumux::start<TypeTag>(argc, argv, usage);
#else
#warning ALUGrid needed for this test.
    std::cout << "ALUGrid needed for this test. Aborting." << std::endl;
#endif
}
