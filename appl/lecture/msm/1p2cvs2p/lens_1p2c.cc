/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder, Klaus Mosthaf                  *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#include "config.h"

#include "lensproblem1p2c.hh"

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

////////////////////////
// helper class for grid instantiation
////////////////////////
template <class Grid, class Scalar>
class CreateGrid
{
};

#if HAVE_UG
template <class Scalar>
class CreateGrid<Dune::UGGrid<2>, Scalar>
{
public:
    static Dune::UGGrid<2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                   const Dune::FieldVector<int, 2> &cellRes)
    {
        Dune::UGGrid<2> *grid = new Dune::UGGrid<2>;
        Dune::GridFactory<Dune::UGGrid<2> > factory(grid);
        for (int i=0; i<=cellRes[0]; i++) {
            for (int j=0; j<=cellRes[1]; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = upperRight[0]*double(i)/cellRes[0];
                pos[1] = upperRight[1]*double(j)/cellRes[1];
                factory.insertVertex(pos);
            }
        }

        for (int i=0; i<cellRes[0]; i++) {
            for (int j=0; j<cellRes[1]; j++) {
                std::vector<unsigned int> v(4);
                v[0] = i*(cellRes[1]+1) + j;
                v[1] = i*(cellRes[1]+1) + j+1;
                v[2] = (i+1)*(cellRes[1]+1) + j;
                v[3] = (i+1)*(cellRes[1]+1) + j+1;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
            }
        }

        factory.createGrid();
        grid->loadBalance();
        return grid;
    }
};
#endif

template <class Scalar>
class CreateGrid<Dune::YaspGrid<2>, Scalar>
{
public:
    static Dune::YaspGrid<2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                     const Dune::FieldVector<int, 2> &cellRes)
    {
        return new Dune::YaspGrid<2>(
#ifdef HAVE_MPI
            Dune::MPIHelper::getCommunicator(),
#endif
            upperRight, // upper right
            cellRes, // number of cells
            Dune::FieldVector<bool,2>(false), // periodic
            4); // overlap
    };
};

template <class Scalar>
class CreateGrid<Dune::SGrid<2, 2>, Scalar>
{
public:
    static Dune::SGrid<2, 2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                     const Dune::FieldVector<int, 2> &cellRes)
    {
        return new Dune::SGrid<2,2>(cellRes, // number of cells
                                    Dune::FieldVector<Scalar, 2>(0.0), // lower left
                                    upperRight); // upper right

    };
};

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s InputFileName\n")%progname;
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        typedef TTAG(LensProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc != 2)
            usage(argv[0]);

        std::string inputFileName;
        inputFileName = argv[1];

        Dune::ParameterTreeParser::readINITree(inputFileName, Params::tree());

        double tEnd, dt;
        tEnd = Params::tree().get<double>("tEnd");
        dt = Params::tree().get<double>("dt");

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        GlobalPosition lowerLeft(0.0);
        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 5.0;
        upperRight[1] = 4.0;
        res[0]        = 40;
        res[1]        = 64;

        std::auto_ptr<Grid> grid(CreateGrid<Grid, Scalar>::create(upperRight, res));

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        // specify dimensions of the low-permeable lens
        GlobalPosition lowerLeftLens, upperRightLens;
        lowerLeftLens[0] = 0.8;
        lowerLeftLens[1] = 2.0;
        upperRightLens[0] = 4.0;
        upperRightLens[1] = 3.0;

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, grid->leafView(), lowerLeft, upperRight, lowerLeftLens, upperRightLens);
        timeManager.init(problem, 0, dt, tEnd);
//        if (restart)
//            problem.restart(restartTime);
        timeManager.run();
        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
