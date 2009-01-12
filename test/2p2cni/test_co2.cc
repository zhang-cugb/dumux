#include "config.h"
#include <iostream>
#define DUMMY
//#undef DUMMY
#ifdef DUMMY
//#ifdef HAVE_UG
#include <iomanip>
//#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "brineco2problem.hh"
//#include "benchmarkproblem.hh"
#include "co2_soilproperties.hh"

#include "dumux/material/phaseproperties/phaseproperties_brineco2.hh"
#include "dumux/material/matrixproperties.hh"
#include "dumux/material/twophaserelations.hh"

#include "dumux/2p2cni/fv/boxco2.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/io/readstarformat.cc"
#include "dumux/io/vtkmultiwriter.hh"


int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;
    double depthBOR = 800.0;  // bottom of reservoir
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(60);
    outerUpperRight[1] = 40;

    if (argc != 4) {
      std::cout << "usage: 2p2cni basefilename tEnd dt" << std::endl;
      return 0;
    }
    std::string arg1(argv[2]);
    std::istringstream is1(arg1);
    double tEnd;
    is1 >> tEnd;
    std::string arg2(argv[3]);
    std::istringstream is2(arg2);
    double dt;
    is2 >> dt;

    // create a grid object
    typedef Dune::SGrid<dim,dim> GridType;
    //typedef Dune::YaspGrid<dim,dim> GridType;
   //  typedef Dune::UGGrid<dim> GridType;
//    typedef Dune::ALUSimplexGrid<dim,dim> GridType;
//    typedef Dune::ALUCubeGrid<dim,dim> GridType;

    Dune::GridPtr<GridType> gridPointer(argv[1]);
    GridType& grid = *gridPointer;
    //readStarFormat(grid, argv[1]);
    //grid.createLGMGrid(argv[1]);

     Dune::gridinfo(grid);

     // choose fluids and properties
     Dune::Liq_BrineCO2 wPhase;
     Dune::Gas_BrineCO2 nPhase;
     Dune::CO2Soil<GridType, NumberType> soil;

     Dune::TwoPhaseRelations<GridType, NumberType>
     materialLaw(soil, wPhase, nPhase);

     Dune::CBrineCO2 multicomp(wPhase, nPhase);

     Dune::BrineCO2Problem<GridType, NumberType> problem(wPhase, nPhase, soil,
                 materialLaw, multicomp, depthBOR);

 //    Dune::BrineCO2Problem<GridType, NumberType> problem(wPhase, nPhase, soil,
//             materialLaw, multicomp, depthBOR);

     typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
     typedef Dune::BoxCO2<GridType, NumberType, MultiWriter> TwoPhase;
      TwoPhase twoPhase(grid, problem);

      Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt, "co2", 1);

      Dune::Timer timer;
      timer.reset();
      MultiWriter writer("co2");

      //timeloop.execute(twoPhase);

//  for timeloop.executeMultiWriter(twoPhase, writer, true) initial
//  values are read from restart file data.dgf
//  at the moment this only works for SGrid in 2D and for ALUCubeGrid in 3D
      timeloop.executeMultiWriter(twoPhase, writer);
      std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

      //std::cout << twoPhase.injected() << " kg CO2 injected." << std::endl;

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
#else

int main (int argc , char **argv) try
{
  //std::cout << "Please install the UG library." << std::endl;
  std::cout << "Dummy implementation, this test would not compile at the moment." << std::endl;

  return 1;
}
catch (...)
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif
