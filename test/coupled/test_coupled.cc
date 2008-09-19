#include "config.h"
#include <iostream>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
//#include <../../../dune-subgrid/subgrid/subgrid.hh>
#include <dune/disc/stokes/dgstokes.hh>

#include "dumux/operators/p1operatorextended.hh"

#include <dumux/timedisc/timeloop.hh>
#include <dumux/coupled/coupledmodel.hh>
#include "boxdiffusion.hh"

#define DGF

namespace Dune
{
template<int dim>
struct P1Layout
{
	bool contains (Dune::GeometryType gt)
	{
		return gt.dim() == 0;
	}
}; 
template<class Grid, class Solution, class Problem> 
double discreteError(const Grid& grid, const Solution& solution, const Problem& problem)
{
	  enum{dim=Grid::dimension};
		typedef typename Grid::LeafGridView GV;
	    typedef typename GV::IndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VM;
		typedef typename GV::template Codim<dim>::Iterator VertexIterator;
	  
	  VM vertexMapper(grid, grid.leafIndexSet());
	  double error = 0.0;
	  const GV& gridview(grid.leafView());
	  
	  VertexIterator endIt = gridview.template end<dim>();
	  VertexIterator it = gridview.template begin<dim>();
	  for (; it != endIt; ++it)
	  {
		  // get exact solution at vertex
		  FieldVector<double,dim> globalCoord = (*it).geometry()[0];
		  double exact = problem.exact(globalCoord);

		  // get approximate solution at vertex
		  int globalId = vertexMapper.map(*it);
		  double approximate = (*solution)[globalId];
		  
		  error += (exact - approximate)*(exact - approximate);
	  }
		  
	  return sqrt(error);
}
}

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    if (argc != 2 && argc != 3) {
      std::cout << "usage: test_coupled dgffilename/basefilename [refinementsteps]" << std::endl;
      return 1;
    }
    int refinementSteps = 0;
    if (argc == 3) {
    	std::string arg2(argv[2]);
    	std::istringstream is2(arg2);
    	is2 >> refinementSteps;
    }
    
    typedef Dune::SGrid<dim,dim> GridType; 
    Dune::GridPtr<GridType> gridPtr( argv[1] );
    GridType& grid = *gridPtr;
    
    if (refinementSteps)
    	grid.globalRefine(refinementSteps);

    Dune::gridinfo(grid);

//    Dune::SubGrid subGrid(grid);
    
    DiffusionParameters<GridType,NumberType> problem;
    
    typedef Dune::LeafP1BoxDiffusion<GridType, NumberType> Diffusion;
    Diffusion diffusion(grid, problem);
    
    const int vOrder = 2; 
    const int pOrder = 1; 
    DGStokesParameters parameters; 
    Example<dim, NumberType> exactSolution;
    DirichletBoundary<GridType> dirichletBoundary(exactSolution);
    RightHandSide<GridType> rightHandSide(exactSolution);
    typedef Dune::DGStokes<GridType, vOrder, pOrder> DGStokes;
	DGStokes dGStokes(grid, exactSolution, parameters, dirichletBoundary, rightHandSide, refinementSteps); 

    typedef Dune::CoupledModel<Diffusion,DGStokes> CoupledModel;
    CoupledModel coupledModel(grid, diffusion, grid, dGStokes, true);
    
    coupledModel.initial();
    coupledModel.assemble();
    printmatrix(std::cout, *(diffusion.A), "local stiffness matrix", "row", 11, 4);
    printmatrix(std::cout, coupledModel.matrix(), "global stiffness matrix", "row", 11, 4);
    
    Dune::TimeLoop<GridType, Diffusion> timeloop(0, 1, 1, "box3d", 1);
    
    timeloop.execute(diffusion);

    printmatrix(std::cout, *(diffusion.A), "local stiffness matrix", "row", 11, 4);
	std::cout << "discrete error = " << discreteError(grid, *diffusion, problem) << std::endl;
	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
