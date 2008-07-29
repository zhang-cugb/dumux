template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
class PipeFlow {
public:
	BCP boundaryConditionP; //!< functor for the boundary condition
	BCV boundaryConditionV; //!< functor for the boundary condition
	ICP initialConditionP;  //!< functor for the initial condition pressure
	ICV initialConditionV;  //!< functor for the initial condition velocity
	SST sourceSink;			//!< functor for source/sink term
	Lmbd frictionCoef;		//!< functor for pipe friction coefficient
	double alphaEX;			//!< mass exchange factor for coupling
	double Temp;			//!< gravity vector
	double density;			//!< fluid property
	double viscosity;		//!< fluid property	
	double diameter;		//!< pipe diameter 
	double roughness;		//!< pipe sand grain roughness
	Press pressurePorous;	//!< Pressure of porous media
	Press pressure;   		//!< vector for the values of the pressure 
	Press velocity;      	//!< block vector for the total velocity
	Grid *grid;           	//!< pointer to grid 
	MapperC0 *mapperC0;     //!< pointer to mapper
	MapperNodeGlobalIDtoOnOutIndexType *mapperGlobalNodeIdtoOnPipeNodeIndex;
	MapperNodeGlobalIDtoOnOutIndexType *mapperGlobalNodeIdtoOutPipeNodeIndex;
	VertexVectorOnLineType vertexVectorOnLine;
	VertexVectorOutLineType vertexVectorOutLine;
	
	void PrintVertexVector ();
	//! sets the vector of unknowns to initial values. 
	void SetInitialSolution ();
	// sets Boundary Condition for velocity, pressure boundary condition does not ly on the pressure index
	void SetDrichletBoundary (double t);

	template<class Matrix>
	void MassEquation (unsigned& k, Matrix& A, Press& f, double t, double dt);

	//! \brief Calculate the pressure.
	void IterationStep_Mass (double t, double dt, Press& pressureIt);

	void Iteration (double t, double dt, double maxdef);

	//! calculate the update vector and estimate the time step size.

	void timeLoop(void);

	//! perform one step of a Runge Kutta scheme

	/*! \brief Constructor. 
	 *
	 *  \param size number of degrees of freedom   
	 *  \param g pointer to grid 
	 *  \param m pointer to mapper 
	 *
	 *  Initializes the vectors \a pressure and \a velocity, the boundary condition \a boundaryCondition, 
	 *  and the initial condition \a initialCondition. 
	 *  Sets the pointers \a grid and \a mapper to \a g and \a m. 
	 */
	PipeFlow(int nNodes, Grid *g, MapperC0 *mC0, MapperNodeGlobalIDtoOnOutIndexType *mapGlobalNodeIDtoPipeNodeOnlineIndex, MapperNodeGlobalIDtoOnOutIndexType *mapGlobalNodeIDtoPipeNodeOutlineIndex, VertexVectorOnLineType *vertexVectorOnL, VertexVectorOutLineType *vertexVectorOutL,
			Press pressureP, double alpEx, double tmp, double dens, double kinematicViscosity, double roughn, double diam) 
	: boundaryConditionP(), boundaryConditionV(),  
	initialConditionP(), initialConditionV(), sourceSink(),
	pressure(nNodes), velocity(nNodes)
	{
		vertexVectorOnLine = *vertexVectorOnL;
		vertexVectorOutLine = *vertexVectorOutL;
		grid = g;
		mapperC0 = mC0;
		mapperGlobalNodeIdtoOnPipeNodeIndex = mapGlobalNodeIDtoPipeNodeOnlineIndex;
		mapperGlobalNodeIdtoOutPipeNodeIndex = mapGlobalNodeIDtoPipeNodeOutlineIndex; 
		
		pressurePorous = pressureP;
		alphaEX = alpEx;
		Temp = tmp;
		density = dens;
		viscosity = kinematicViscosity;
		roughness = roughn;
		diameter = diam;
	}
};

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::PrintVertexVector()
{
	
    std::cout << "number of vertices on line = "<< vertexVectorOnLine.size() << std::endl;
    std::cout << "number of vertices out line = "<< vertexVectorOutLine.size() << std::endl;
    
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
    	std::cout << "vertice on line coord: " <<vertexVectorOnLine[k].nodePoint->geometry()[0] << std::endl;
    	std::cout << "       vertice on line: " << k << std::endl;
    	std::cout << "       globalId: " << vertexVectorOnLine[k].globalId << std::endl;
    	std::cout << "       vertexOnLineIndex: " << ((*mapperGlobalNodeIdtoOnPipeNodeIndex).find(vertexVectorOnLine[k].globalId))->second << std::endl;
    	for (unsigned n = 0; n < vertexVectorOnLine[k].parameter.size(); n++)
    	{
    		std::cout << "       parameters: " << vertexVectorOnLine[k].parameter[n] << std::endl;
    	}
    	std::cout << "           boundary: " << vertexVectorOnLine[k].boundary() << std::endl;
    	if (vertexVectorOnLine[k].boundary())
    	{
    		std::cout << "           normalBF: " << vertexVectorOnLine[k].normalBF(vertexVectorOnLine) << std::endl;
    		std::cout << "           unitPDBF: " << vertexVectorOnLine[k].unitPDBF(vertexVectorOnLine) << std::endl;
    	}
    	std::cout << "             length: " << vertexVectorOnLine[k].length(vertexVectorOnLine) << std::endl;
    	std::cout << "             on line size: " << vertexVectorOnLine[k].lineVectorOnLine.size() << std::endl;
        for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOnLine.size(); m++)
        {
        	std::cout << "                neighbour vertices on line: " << vertexVectorOnLine[k].indexVertexVectorOnLine[m] << std::endl;
        	std::cout << "                normal vector: " << vertexVectorOnLine[k].normal(vertexVectorOnLine, vertexVectorOnLine[k].indexVertexVectorOnLine[m]) << std::endl;
        	std::cout << "                unitPD vector: " << vertexVectorOnLine[k].unitPD(vertexVectorOnLine, vertexVectorOnLine[k].indexVertexVectorOnLine[m]) << std::endl;
        }
    	std::cout << "             out line size: " << vertexVectorOnLine[k].lineVectorOutLine.size() << std::endl;
        for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOutLine.size(); m++)
        {
        	std::cout << "                neighbour vertices out line: " << vertexVectorOnLine[k].indexVertexVectorOutLine[m] << std::endl;
        	std::cout << "                normal vector: " << vertexVectorOnLine[k].normal(vertexVectorOutLine, vertexVectorOnLine[k].indexVertexVectorOutLine[m]) << std::endl;
        }
    }
    
    for (unsigned k = 0; k < vertexVectorOutLine.size(); k++)
    {
    	std::cout << "vertice out line coord: " <<vertexVectorOutLine[k].nodePoint->geometry()[0] << std::endl;
    	std::cout << "       vertice out line: " << k << std::endl;
    	std::cout << "       globalId: " << vertexVectorOutLine[k].globalId << std::endl;
    	std::cout << "       vertexOutLineIndex: " << ((*mapperGlobalNodeIdtoOutPipeNodeIndex).find(vertexVectorOutLine[k].globalId))->second << std::endl;
    	std::cout << "          out line size: " << vertexVectorOutLine[k].lineVectorOutLine.size() << std::endl;
        for (unsigned m = 0; m < vertexVectorOutLine[k].indexVertexVectorOnLine.size(); m++)
        {
        	std::cout << "             neighbour vertices on line: " << vertexVectorOutLine[k].indexVertexVectorOnLine[m] << std::endl;
        }
    }

}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::SetInitialSolution ()
{
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
    	typedef Dune::FieldVector<double,Grid::dimension>  FieldVector;
    	FieldVector global = vertexVectorOnLine[k].nodePoint->geometry()[0];
		pressure[k] = initialConditionP(global);
    }
	return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::SetDrichletBoundary (double t)
{
	
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
    	if (vertexVectorOnLine[k].boundary() )
    	{
    		int boundaryId = (int) vertexVectorOnLine[k].parameter[1]; 
        	int boundaryType;
           	double boundaryValue = boundaryConditionP(boundaryId, t, boundaryType);
           	if (boundaryType == 1) // if dirichlet boundary
           	{
               	pressure[k]=boundaryValue;
           	}
    	}
    }
	return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
template<class Matrix>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::MassEquation(unsigned& k, Matrix& A, Press& f, double t, double dt)
{
	  // first we extract the dimensions of the grid  
	  const int dim = Grid::dimension;

	  int indexi = k;
	  
	  double length =  vertexVectorOnLine[k].length(vertexVectorOnLine);
	  double PI = 3.14;	  
	  double crossArea = PI* diameter * diameter / 4;
	  
	  double TaoNoVel= 8 * viscosity * density /diameter * (PI * diameter) / crossArea; // Tao(no velocity inside) * perimeter / crossarea uniform
	  double K = 1 / TaoNoVel; // at the face
	  
	  typedef Dune::FieldVector<double,dim> FieldVector;
	  FieldVector globali = vertexVectorOnLine[k].nodePoint->geometry()[0];
	  
	  f[indexi] += sourceSink(globali, t) / density *length*crossArea;

	  f[indexi] += alphaEX * pressurePorous[vertexVectorOnLine[k].globalId];
	  

	  double prepor = pressurePorous[vertexVectorOnLine[k].globalId];
	
	  std::cout<<"globalID "<< vertexVectorOnLine[k].globalId << " pressurePorous " << prepor <<std::endl;
	  
	  A[indexi][indexi] += alphaEX * 1.0;
	  
	  for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOnLine.size(); m++)
	  {
		  int indexj = vertexVectorOnLine[k].indexVertexVectorOnLine[m];

		  FieldVector globalj = vertexVectorOnLine[indexj].nodePoint->geometry()[0];
		  FieldVector distVect = globalj;
		  distVect -= globali;
		  double distVal = distVect.two_norm();

		  typedef Dune::FieldVector<double,dim> FieldVector;
		  FieldVector unitOuterNormal = vertexVectorOnLine[k].normal(vertexVectorOnLine, indexj);
		  
		  double sign = distVect * unitOuterNormal;
		  if (sign>0) sign=1.0;
		  else sign=-1.0; 
		  
		  A[indexi][indexi] += 1.0 * K * crossArea * 1.0 * sign / distVal;
		  A[indexi][indexj] += -1.0 * K * crossArea * 1.0 * sign / distVal;
	  }

	  // treat Neumann boundary  
	  if (vertexVectorOnLine[k].boundary() )
	  {
		  int boundaryId = (int) vertexVectorOnLine[k].parameter[1]; 
		  int boundaryType;
		  double boundaryValue = boundaryConditionP(boundaryId, t, boundaryType);
		  if (boundaryType == 2) // if neumann flow boundary
		  {

			  f[indexi] += -1 * crossArea * boundaryValue; // positiv boundaryValue --> out flow,  negativ boundaryValue --> in flow
		  }
	  }

	return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::timeLoop(void)
{
	double tstart;
	double tend;
	double max_dt; 
	double first_dt; 
	double CFL_factor;
	int flag;
	int n_iter;
	double max_def;
	int modulo, stages;

	// get timeloop options
	TimeloopOptsPipe( tstart, tend, max_dt, first_dt, CFL_factor, flag, n_iter, max_def, modulo, stages );

	// initialize concentration with initial values
	SetInitialSolution();  

	// printvector(std::cout,pressure,"initial pressure","row",200,4);

	//  vtkout(*grid, pressure, permeability, velocity, "multiscalec", 0);
	//printvector(std::cout,velocity,"velocity","row",4,4);  

	// generate one meta vtk-file holding the individual timesteps 
	char fileName[128] = "pipe";
	char multiFileName[128];	
	char fileNameVTK[128];	
	sprintf(multiFileName,"multi-%s.pvd",fileName );
	std::ofstream multiFile(multiFileName);
	multiFile << "<?xml version=\"1.0\"?>" << std::endl 
	<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" " << std::endl
	<< "compressor=\"vtkZLibDataCompressor\">"  << std::endl
	<< " <Collection>" << std::endl;
	
	// now do the time steps
	double dt = max_dt; 
	double t = tstart;
	int k = 0;
	std::cout.setf (std::ios::scientific, std::ios::floatfield);
	while (t<tend+dt/2.0) {
		
		SetDrichletBoundary (t);

		// outputs
		std::cout << "timestep: " << k << "\t t=" << t << "\t dt=" << dt << std::endl;
		std::cout.setf (std::ios::scientific, std::ios::floatfield);
		// printvector(std::cout,pressure,"pressure","row",200,1);
		
//		printvector(std::cout,pressure,"pressure","row",200,1);
//		printvector(std::cout,velocity,"velocity","row",200,1);

		if (k%modulo==0) 
		{
//			vtkout_pipeflow(*grid, pressure, velocity, fileName,k/modulo); // TODO
			
			sprintf(fileNameVTK,"%s-%05d.vtp", fileName, k/modulo);
			multiFile << "   <DataSet timestep=\"" << t << "\" file=\"" 
			<< fileNameVTK << "\"/>" << std::endl;			 
		}

		{ 
			Iteration (t, dt, max_def); 
		}

		k++;
		t += dt;
	} 
	  // finalize output
	  multiFile << " </Collection>" << std::endl << "</VTKFile>" << std::endl;
	  multiFile.close();
	  
	return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::Iteration (double t, double dt, double max_def)
{
//	int nElem =  vertexVectorOnLine.size();
//	int systemSize = nElem;

	Press pressureIt = pressure;
	
	for (int n = 0; n < 101; n++)
	{
		Press pressureOldIt = pressureIt;

		IterationStep_Mass(t, dt, pressureIt);  
//		printvector(std::cout, pressureIt, "PressureIt","row",100,1,4);

		Press deltaPressure = pressureIt;
		deltaPressure -= pressureOldIt;
		
		double deltaPressureTwoNorm = deltaPressure.two_norm();
		double pressureItTwoNorm = pressureIt.two_norm();
		double defectPressure;

		if (pressureItTwoNorm == 0.0)
		{
			defectPressure = 0.0;
		}
		else
		{
			defectPressure = deltaPressureTwoNorm/pressureItTwoNorm;
		}
		
		std::cout << "defectPressure=" << defectPressure << std::endl;

		if (defectPressure < max_def)
		{
			pressure= pressureIt;
//			std::cout << "Converged in  n= " << n << "iterations" << std::endl;
//			printvector(std::cout,pressure,"Pressure","row",100,1,4);
			break;
		}
	    else if (n >= 100)
	    {
	    	Dune::Exception exception;
	    	exception.message("Diverged !!!!");
	    	throw exception;
	    }
	}
	
  return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class Grid, class MapperC0, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, Grid, MapperC0, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::IterationStep_Mass (double t, double dt, Press& pressureIt)
{
	int nElem = vertexVectorOnLine.size();
	int systemSize = nElem;

	// solution vector SolutionVector and right side vector f
	Press f(systemSize);
	Press SolutionVector(systemSize);
	f=0;
	SolutionVector=0;
	typedef Dune::FieldMatrix<double,1,1> MB;
	Dune::BCRSMatrix<MB> A(systemSize, systemSize, Dune::BCRSMatrix<MB>::random);
	
	// determine matrix row sizes 
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
		// cell index
		int indexi = k;
		
		// initialize row size = diagonal 
		int rowSize = systemSize;

		A.setrowsize(indexi, rowSize);
    }
	A.endrowsizes();

	// determine position of matrix entries 
	for (int i = 0;i<systemSize; ++i)
	{
		for (int j = 0;j<systemSize; ++j)
		{
			A.addindex(i, j);
		}
	}
	A.endindices();		

	//initiliaze A with 0
	for (int i = 0;i<systemSize; ++i)
	{	
		f[i]= 0.0;
		for (int j = 0;j<systemSize; ++j)
		{
			A[i][j]= 0.0;
		}
	}

	// fill matrix
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
    	MassEquation<Dune::BCRSMatrix<MB> > (k, A, f, t, dt);
    }   // end grid traversal
    
    // set dirichlet boundary in Matrix A and right hand side f
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
    {
    	if (vertexVectorOnLine[k].boundary() )
    	{
    		int boundaryId = (int) vertexVectorOnLine[k].parameter[1]; 
        	int boundaryType;
           	double boundaryConditionValue = boundaryConditionP(boundaryId, t, boundaryType);
           	if (boundaryType == 1) // if dirichlet boundary
           	{
           		int indexi = k;
           		A[indexi]=0;
           		A[indexi][indexi]=1.0;
           		f[indexi]= boundaryConditionValue;
           	}
    	}
    }   // end grid traversal    

//	printvector(std::cout,f,"right hand side","row",200,1);
//	printmatrix(std::cout,A,"matrix","",8,1);
	// set up the high-level solver objects
	typedef Dune::FieldVector<double, 1> VB;
	typedef Dune::BlockVector<VB> Vector;
	typedef Dune::BCRSMatrix<MB> Matrix;
	Dune::MatrixAdapter<Matrix,Vector,Vector> op(A);        // make linear operator from A
	Dune::InverseOperatorResult r;
	
    Dune::SeqPardiso<Matrix,Vector,Vector> pardiso(A);        // preconditioner object
    Dune::LoopSolver<Vector> loop(op, pardiso, 1E-14, 2, 1);  // an inverse operator
	loop.apply(SolutionVector, f, r);  						  // call the solver
	
//	Dune::SeqGS<Matrix,Vector,Vector> ilu0(bigMatrix, 1, 1.0);         // preconditioner object
//	Dune::BiCGSTABSolver<Vector> cg(op, ilu0, 1E-14, 1000, 2); // an inverse operator
//	cg.apply(bigSolutionVector, fBig, r);							   // call the solver	
	

	//fill deltaVelocity 
	for (int i = 0;i<nElem; ++i)
	{					
		pressureIt[i]=SolutionVector[i];
	}

//	printvector(std::cout,deltaVelocity,"deltaVelocity","row",200,1);
//	printvector(std::cout,deltaPressure,"deltaPressure","row",200,1);
	return;
}
