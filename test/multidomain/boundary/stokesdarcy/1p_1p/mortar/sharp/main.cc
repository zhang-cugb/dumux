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
 * \ingroup BoundaryTests
 * \brief TODO doc me
 */
#include <config.h>
#include <iostream>
#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/integrate.hh>
#include <dumux/discretization/normalfluxbasis.hh>
#include <dumux/discretization/fem/fegridgeometry.hh>

#include "../problem_darcy.hh"
#include "../problem_stokes.hh"

#include "../subdomainsolvers.hh"

#include "projector.hh"
#include "../reconstructor.hh"

#include "preconditioner.hh"
#include "interfaceoperator.hh"

////////////////////////////////////////////
// Some aliases etc to be used in solve() //
////////////////////////////////////////////

// Traits class for a sub-domain
template<class Solver>
struct SubDomainTraits
{
    using SolutionVector = typename Solver::SolutionVector;
    using GridGeometry = typename Solver::FVGridGeometry;
    using GridVariables = typename Solver::GridVariables;
    using FluxVariables = typename Solver::FluxVariables;
};

// Define grid and basis for mortar domain
struct MortarTraits
{
    static constexpr int basisOrder = 1;

    using Scalar = double;
    using Grid = Dune::FoamGrid<1, 2>;
    using GridView = typename Grid::LeafGridView;
    using BlockType = Dune::FieldVector<Scalar, 1>;
    using SolutionVector = Dune::BlockVector<BlockType>;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, basisOrder>;
    using GridGeometry = Dumux::FEGridGeometry<FEBasis>;
};

template<class SubDomainTypeTag>
using DarcySolverType = Dumux::DarcySolver<SubDomainTypeTag>;

template<class SubDomainTypeTag>
using StokesSolverType = Dumux::StokesSolver<SubDomainTypeTag>;

/////////////////////////////////
// The iterative solve routine //
/////////////////////////////////
template< class Solver1,
          class Solver2 >
void solveMortar()
{
    using namespace Dumux;

    Dune::Timer watch;

    // create sub-domain solvers
    auto solver1 = std::make_shared< Solver1 >();
    auto solver2 = std::make_shared< Solver2 >();

    solver1->init("Domain1");
    solver2->init("Domain2");

    // make mortar grid, function space basis and solution
    using MortarGrid = typename MortarTraits::Grid;
    using MortarGridView = typename MortarTraits::GridView;
    GridManager<MortarGrid> mortarGridManager;
    mortarGridManager.init("Mortar");

    const auto& mortarGridView = mortarGridManager.grid().leafGridView();
    auto feBasis = std::make_shared<typename MortarTraits::FEBasis>(mortarGridView);

    using MortarGridGeometry = typename MortarTraits::GridGeometry;
    auto mortarGridGeometry = std::make_shared<MortarGridGeometry>(feBasis);

    using MortarSolution = typename MortarTraits::SolutionVector;
    auto mortarSolution = std::make_shared<MortarSolution>();
    mortarSolution->resize(feBasis->size());
    *mortarSolution = 0.0;

    // create vtk writer for mortar grid
    auto mortarWriter = std::make_shared<Dune::VTKWriter<MortarGridView>>(mortarGridView, Dune::VTK::nonconforming);
    Dune::VTKSequenceWriter<MortarGridView> mortarSequenceWriter(mortarWriter, "mortar_sharp");

    auto mortarGridFunction = Dune::Functions::template makeDiscreteGlobalBasisFunction<typename MortarSolution::block_type>(*feBasis, *mortarSolution);
    const auto fieldInfoMortar = Dune::VTK::FieldInfo({"flux", Dune::VTK::FieldInfo::Type::scalar, 1});

    if (MortarTraits::basisOrder == 0)
        mortarWriter->addCellData(mortarGridFunction, fieldInfoMortar);
    else
        mortarWriter->addVertexData(mortarGridFunction, fieldInfoMortar);

    // project initial mortar solution into sub-domains
    const auto glue1 = makeGlue(*solver1->gridGeometryPointer(), *mortarGridGeometry);
    const auto glue2 = makeGlue(*solver2->gridGeometryPointer(), *mortarGridGeometry);

    const auto& basis1 = getFunctionSpaceBasis(*solver1->gridGeometryPointer());
    const auto& basis2 = getFunctionSpaceBasis(*solver2->gridGeometryPointer());

    const auto matrices1 = makeProjectionMatricesPair(basis1, getFunctionSpaceBasis(*mortarGridGeometry), glue1).second;
    const auto matrices2 = makeProjectionMatricesPair(basis2, getFunctionSpaceBasis(*mortarGridGeometry), glue2).second;
    const SharpMortarProjector<MortarSolution> projector(matrices1.first, matrices1.second, matrices2.first, matrices2.second);

    auto x1 = *mortarSolution;
    auto x2 = *mortarSolution;

    if (solver1->problemPointer()->isOnNegativeMortarSide()) x1 *= -1.0;
    if (solver2->problemPointer()->isOnNegativeMortarSide()) x2 *= -1.0;

    const auto p = projector.projectMortarToSubDomain(x1, x2);

    using namespace Dune::Indices;
    solver1->problemPointer()->setMortarProjection(p[_0]);
    solver2->problemPointer()->setMortarProjection(p[_1]);

    // write out initial solution
    mortarSequenceWriter.write(0.0);
    solver1->write(0.0);
    solver2->write(0.0);

    // create interface operator
    using Reconstructor1 = MortarReconstructor< SubDomainTraits<Solver1> >;
    using Reconstructor2 = MortarReconstructor< SubDomainTraits<Solver2> >;
    using Operator = SharpMortarInterfaceOperator<MortarSolution, MortarGridGeometry,
                                                  Solver1, Reconstructor1,
                                                  Solver2, Reconstructor2 >;
    Operator op(solver1, solver2, mortarGridGeometry);

    // first compute the jump in mortar variable
    solver1->problemPointer()->setUseHomogeneousSetup(false);
    solver2->problemPointer()->setUseHomogeneousSetup(false);

    MortarSolution deltaMortarVariable;
    op.apply(*mortarSolution, deltaMortarVariable);

    // Solve the homogeneous problem with CG solver
    const double reduction = getParam<double>("InterfaceSolver.ResidualReduction");
    const std::size_t maxIt = getParam<double>("InterfaceSolver.MaxIterations");
    const int verbosity = getParam<int>("InterfaceSolver.Verbosity");

    solver1->problemPointer()->setUseHomogeneousSetup(true);
    solver2->problemPointer()->setUseHomogeneousSetup(true);

    // create preconditioner
    using Prec = Dumux::SharpMortarPreconditioner<MortarSolution, MortarGridGeometry,
                                                  Solver1, Reconstructor1,
                                                  Solver2, Reconstructor2 >;
    Prec prec(solver1, solver2, mortarGridGeometry);

    // apply linear solver using our linear operator
    deltaMortarVariable *= -1.0;
    Dune::InverseOperatorResult result;

    const auto lsType = getParam<std::string>("InterfaceSolver.LinearSolverType");
    if (lsType == "CG")
    {
        Dune::CGSolver<MortarSolution> cgSolver(op, prec, reduction, maxIt, verbosity);
        cgSolver.apply(*mortarSolution, deltaMortarVariable, result);
    }
    else if (lsType == "GMRes")
    {
        Dune::RestartedGMResSolver<MortarSolution> gmresSolver(op, prec, reduction, maxIt, maxIt, verbosity);
        gmresSolver.apply(*mortarSolution, deltaMortarVariable, result);
    }

    if (!result.converged)
        DUNE_THROW(Dune::InvalidStateException, "CG solver did not converge with given maximum number of iterations");

    // solve the sub-domains again to get the right output
    solver1->problemPointer()->setUseHomogeneousSetup(false);
    solver2->problemPointer()->setUseHomogeneousSetup(false);
    op.apply(*mortarSolution, deltaMortarVariable);

    // add the recovered pressures from the sub-domain to the vtk output
    using R1 = Reconstructor1;
    using R2 = Reconstructor2;
    auto p1 = R1::template recoverSolution<MortarSolution>(*solver1->gridGeometryPointer(),
                                                           *solver1->gridVariablesPointer(),
                                                           *solver1->solutionPointer(),
                                                            op.coupledScvfMap1());

    auto p2 = R2::template recoverSolution<MortarSolution>(*solver2->gridGeometryPointer(),
                                                           *solver2->gridVariablesPointer(),
                                                           *solver2->solutionPointer(),
                                                            op.coupledScvfMap2());

    solver1->outputModule().addField(p1, "ifPressure");
    solver2->outputModule().addField(p2, "ifPressure");

    // add interface pressure projected onto the mortar domain
    auto fluxMortar = [&] (const auto& pos) { return solver1->problemPointer()->exactFlux(pos)[1]; };
    auto pressureTuple = op.projector().projectSubDomainToMortar(p1, p2);
    const auto gfProject1 = Dune::Functions::template makeDiscreteGlobalBasisFunction<typename MortarSolution::block_type>(*feBasis, pressureTuple[_0]);
    const auto gfProject2 = Dune::Functions::template makeDiscreteGlobalBasisFunction<typename MortarSolution::block_type>(*feBasis, pressureTuple[_1]);
    const auto analyticFluxMortar = Dune::Functions::makeAnalyticGridViewFunction(fluxMortar, feBasis->gridView());
    const auto fieldInfoProjection1 = Dune::VTK::FieldInfo({"ifPressureProjected1", Dune::VTK::FieldInfo::Type::scalar, 1});
    const auto fieldInfoProjection2 = Dune::VTK::FieldInfo({"ifPressureProjected2", Dune::VTK::FieldInfo::Type::scalar, 1});

    // project mortar projections from sub-domains back to mortar
    const auto projectedFlux1 = solver1->problemPointer()->mortarProjection();
    const auto projectedFlux2 = solver2->problemPointer()->mortarProjection();
    auto fluxTuple = op.projector().projectSubDomainToMortar(projectedFlux1, projectedFlux2);
    const auto gfBackProjectFlux1 = Dune::Functions::template makeDiscreteGlobalBasisFunction<typename MortarSolution::block_type>(*feBasis, fluxTuple[_0]);
    const auto gfBackProjectFlux2 = Dune::Functions::template makeDiscreteGlobalBasisFunction<typename MortarSolution::block_type>(*feBasis, fluxTuple[_1]);
    const auto fieldInfoBackProjection1 = Dune::VTK::FieldInfo({"fluxBackProjection1", Dune::VTK::FieldInfo::Type::scalar, 1});
    const auto fieldInfoBackProjection2 = Dune::VTK::FieldInfo({"fluxBackProjection2", Dune::VTK::FieldInfo::Type::scalar, 1});

    if (MortarTraits::basisOrder == 0)
    {
        mortarWriter->addCellData(gfProject1, fieldInfoProjection1);
        mortarWriter->addCellData(gfProject2, fieldInfoProjection2);
        mortarWriter->addCellData(gfBackProjectFlux1, fieldInfoBackProjection1);
        mortarWriter->addCellData(gfBackProjectFlux2, fieldInfoBackProjection2);
        mortarWriter->addCellData(analyticFluxMortar, Dune::VTK::FieldInfo("flux_exact", Dune::VTK::FieldInfo::Type::scalar, 1));
    }
    else
    {
        mortarWriter->addVertexData(gfProject1, fieldInfoProjection1);
        mortarWriter->addVertexData(gfProject2, fieldInfoProjection2);
        mortarWriter->addVertexData(gfBackProjectFlux1, fieldInfoBackProjection1);
        mortarWriter->addVertexData(gfBackProjectFlux2, fieldInfoBackProjection2);
        mortarWriter->addVertexData(analyticFluxMortar, Dune::VTK::FieldInfo("flux_exact", Dune::VTK::FieldInfo::Type::scalar, 1));
    }

    // write solutions
    mortarSequenceWriter.write(1.0);
    solver1->write(1.0);
    solver2->write(1.0);

    // compute pressure L2 error norm
    using namespace Dune::Functions;
    using BlockType = typename MortarSolution::block_type;
    auto gf1 = makeDiscreteGlobalBasisFunction<BlockType>(basis1, *solver1->solutionPointer());
    auto gf2 = makeDiscreteGlobalBasisFunction<BlockType>(basis2, *solver2->solutionPointer());

    auto pressure1 = [&] (const auto& pos) { return solver1->problemPointer()->exact(pos); };
    auto pressure2 = [&] (const auto& pos) { return solver2->problemPointer()->exact(pos); };
    auto analytic1 = makeAnalyticGridViewFunction(pressure1, basis1.gridView());
    auto analytic2 = makeAnalyticGridViewFunction(pressure2, basis2.gridView());

    const int order = getParam<int>("L2Error.IntegrationOrder");
    const auto l2norm1 = Dumux::integrateL2Error(basis1.gridView(), analytic1, gf1, order);
    const auto l2norm2 = Dumux::integrateL2Error(basis2.gridView(), analytic2, gf2, order);
    std::cout << "Pressure norms: " << l2norm1 << " - " << l2norm2 << std::endl;

    // compute flux L2 error norm
    const auto fluxBasis1 = getVelocityFunctionSpaceBasis(*solver1->gridGeometryPointer());
    const auto fluxBasis2 = getVelocityFunctionSpaceBasis(*solver2->gridGeometryPointer());

    const auto coeff1 = getVelocityCoefficientVector<typename Solver1::FluxVariables>(fluxBasis1, *solver1->gridGeometryPointer(), *solver1->gridVariablesPointer(), *solver1->solutionPointer());
    const auto coeff2 = getVelocityCoefficientVector<typename Solver2::FluxVariables>(fluxBasis2, *solver2->gridGeometryPointer(), *solver2->gridVariablesPointer(), *solver2->solutionPointer());

    using FluxRange = Dune::FieldVector<typename MortarTraits::Scalar, int(MortarGridView::dimension) + 1>;
    auto gfFlux1 = makeDiscreteGlobalBasisFunction<FluxRange>(fluxBasis1, coeff1[/*phaseIdx*/0]);
    auto gfFlux2 = makeDiscreteGlobalBasisFunction<FluxRange>(fluxBasis2, coeff2[/*phaseIdx*/0]);
    auto gfFluxMortar = makeDiscreteGlobalBasisFunction<typename MortarTraits::BlockType>(*feBasis, *mortarSolution);

    auto flux1 = [&] (const auto& pos) { return solver1->problemPointer()->exactFlux(pos); };
    auto flux2 = [&] (const auto& pos) { return solver2->problemPointer()->exactFlux(pos); };
    auto analyticFlux1 = makeAnalyticGridViewFunction(flux1, fluxBasis1.gridView());
    auto analyticFlux2 = makeAnalyticGridViewFunction(flux2, fluxBasis2.gridView());

    Dune::VTKWriter<typename Solver1::FVGridGeometry::GridView> fluxWriter1(solver1->gridGeometryPointer()->gridView());
    Dune::VTKWriter<typename Solver2::FVGridGeometry::GridView> fluxWriter2(solver2->gridGeometryPointer()->gridView());

    // write out exact interface flux as cell data with entries only in the interface cells
    MortarSolution exactIFFlux1; exactIFFlux1.resize(solver1->gridGeometryPointer()->gridView().size(0));
    MortarSolution exactIFFlux2; exactIFFlux2.resize(solver2->gridGeometryPointer()->gridView().size(0));

    for (const auto& e : elements(solver1->gridGeometryPointer()->gridView()))
    {
        auto fvGeometry = localView(*solver1->gridGeometryPointer());
        fvGeometry.bind(e);

        for (const auto& scvf : scvfs(fvGeometry))
            if (solver1->problemPointer()->isOnMortarInterface(scvf.ipGlobal()))
                exactIFFlux1[scvf.insideScvIdx()] = solver1->problemPointer()->exactFlux(scvf.ipGlobal())[1];
    }

    for (const auto& e : elements(solver2->gridGeometryPointer()->gridView()))
    {
        auto fvGeometry = localView(*solver2->gridGeometryPointer());
        fvGeometry.bind(e);

        for (const auto& scvf : scvfs(fvGeometry))
            if (solver2->problemPointer()->isOnMortarInterface(scvf.ipGlobal()))
                exactIFFlux2[scvf.insideScvIdx()] = solver2->problemPointer()->exactFlux(scvf.ipGlobal())[1];
    }

    if (solver1->problemPointer()->isOnNegativeMortarSide()) exactIFFlux1 *= -1.0;
    if (solver2->problemPointer()->isOnNegativeMortarSide()) exactIFFlux2 *= -1.0;

    fluxWriter1.addCellData(exactIFFlux1, "exact_if_flux");
    fluxWriter1.addCellData(gfFlux1, Dune::VTK::FieldInfo("flux", Dune::VTK::FieldInfo::Type::vector, 2));
    fluxWriter1.addCellData(analyticFlux1, Dune::VTK::FieldInfo("flux_exact", Dune::VTK::FieldInfo::Type::vector, 2));

    fluxWriter2.addCellData(exactIFFlux2, "exact_if_flux");
    fluxWriter2.addCellData(gfFlux2, Dune::VTK::FieldInfo("flux", Dune::VTK::FieldInfo::Type::vector, 2));
    fluxWriter2.addCellData(analyticFlux2, Dune::VTK::FieldInfo("flux_exact", Dune::VTK::FieldInfo::Type::vector, 2));

    fluxWriter1.write("flux1_sharp");
    fluxWriter2.write("flux2_sharp");

    const auto fluxNnorm1 = Dumux::integrateFaceFluxError(basis1.gridView(), analyticFlux1, gfFlux1, order);
    const auto fluxNnorm2 = Dumux::integrateFaceFluxError(basis2.gridView(), analyticFlux2, gfFlux2, order);
    const auto fluxNormMortar = Dumux::integrateL2Error(feBasis->gridView(), analyticFluxMortar, gfFluxMortar, order);
    std::cout << "Flux norms: " <<  fluxNnorm1 << " - " << fluxNnorm2 << std::endl;

    // lambda to add error from sub-domain
    auto computeIFError = [&] (const auto& solver)
    {
        double fluxProjectionError = 0.0;
        for (const auto& element : elements(solver->gridGeometryPointer()->gridView()))
        {
            auto fvGeometry  = localView(*solver->gridGeometryPointer());
            fvGeometry.bindElement(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (solver->problemPointer()->isOnMortarInterface(scvf.ipGlobal()))
                {
                    const auto discreteFlux = solver->problemPointer()->isOnNegativeMortarSide()
                                              ? -1.0*solver->problemPointer()->mortarProjection()[scvf.insideScvIdx()]
                                              : solver->problemPointer()->mortarProjection()[scvf.insideScvIdx()];
                    const auto geometry = scvf.geometry();
                    for (const auto& ip : Dune::QuadratureRules<double, 1>::rule(geometry.type(), order))
                    {
                        const auto globalPos = geometry.global(ip.position());
                        const auto exactFlux = solver1->problemPointer()->exactFlux(globalPos)*scvf.unitOuterNormal();
                        fluxProjectionError += (discreteFlux-exactFlux)
                                               *(discreteFlux-exactFlux)
                                               *ip.weight()
                                               *geometry.integrationElement(ip.position());
                    }
                }
            }
        }

        using std::sqrt;
        return sqrt(fluxProjectionError);
    };

    const auto ifError1 = computeIFError(solver1);
    const auto ifError2 = computeIFError(solver2);

    // write into file
    using std::sqrt;
    std::ofstream errorFile(getParam<std::string>("L2Error.OutputFile"), std::ios::app);
    errorFile << sqrt(l2norm1*l2norm1+l2norm2*l2norm2) << ","
              << sqrt(fluxNnorm1*fluxNnorm1+fluxNnorm2*fluxNnorm2) << ","
              << fluxNormMortar << ","
              << sqrt(ifError1*ifError1 + ifError2*ifError2) << std::endl;
    errorFile.close();

    // print time necessary for solve
    std::cout << "\n#####################################################\n\n"
              << "Iterative scheme took " << watch.elapsed() << " seconds\n"
              << "\n#####################################################\n\n";
}

///////////////////////////////////////////////////////////////////
// Main Program. Selects the solvers etc to be passed to solve() //
///////////////////////////////////////////////////////////////////
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    if (getParam<std::string>("Mortar.VariableType") != "Flux")
        DUNE_THROW(Dune::InvalidStateException, "Sharp projection scheme only works for flux mortars");

    // get solver types of the two subdomains
    const auto solver1Type = getParam<std::string>("Domain1.SolverType");
    const auto solver2Type = getParam<std::string>("Domain2.SolverType");

    // discretization scheme used in the sub-domains
    const auto discScheme1 = getParam<std::string>("Domain1.DiscretizationScheme");
    const auto discScheme2 = getParam<std::string>("Domain2.DiscretizationScheme");

    //////////////////////////////////////////
    // Check validity of the specifications //
    //////////////////////////////////////////
    for (unsigned int i = 0; i < 2; ++i)
    {
        const auto& s = i == 0 ? solver1Type : solver2Type;
        const auto& d = i == 0 ? discScheme1 : discScheme2;

        if (s != "Darcy") DUNE_THROW(Dune::InvalidStateException, "Sharp projection only works for Darcy-Darcy");
        if (d != "Tpfa") DUNE_THROW(Dune::InvalidStateException, "Sharp projection only works for TPFA");
    }


    // run the main program
    using TTDarcyTpfa = Properties::TTag::DarcyOnePTpfa;
    solveMortar< DarcySolverType<TTDarcyTpfa>, DarcySolverType<TTDarcyTpfa> >();

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
