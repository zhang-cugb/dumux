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
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <iostream>
#include <utility>

#include <dune/common/test/iteratortest.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/discretization/staggered/globalfvgeometry.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/discretization/staggered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/subcontrolvolumeface.hh>

namespace Dumux
{

template<class TypeTag>
class MockProblem
{
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, DofMapper);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    MockProblem(const GridView& gridView) : mapper_(gridView) {}

    const ElementMapper& elementMapper() const
    { return mapper_; }
private:
    ElementMapper mapper_;
};

namespace Properties
{
NEW_TYPE_TAG(TestFVGeometry, INHERITS_FROM(StaggeredModel));

SET_TYPE_PROP(TestFVGeometry, Grid, Dune::YaspGrid<2>);

SET_TYPE_PROP(TestFVGeometry, Problem, Dumux::MockProblem<TypeTag>);

SET_BOOL_PROP(TestFVGeometry, EnableGlobalFVGeometryCache, true);
}

}

template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& t){}
};

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    // aliases
    using TypeTag = TTAG(TestFVGeometry);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename Grid::LeafGridView;

    constexpr int dim = GridView::dimension;
    constexpr int dimworld = GridView::dimensionworld;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GlobalFVGeometry = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{2, 2}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
    auto leafGridView = grid->leafGridView();

    Problem problem(leafGridView);

    GlobalFVGeometry global(leafGridView);
    global.update(problem);

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = problem.elementMapper().index(element);
        std::cout << std::endl << "Checking fvGeometry of element " << eIdx << std::endl;
        auto fvGeometry = localView(global);
        fvGeometry.bind(element);

        auto range = scvs(fvGeometry);
        NoopFunctor<SubControlVolume> op;
        if(0 != testForwardIterator(range.begin(), range.end(), op))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scv : scvs(fvGeometry))
        {
            std::cout << "-- scv " << scv.index() << " center at: " << scv.center() << std::endl;
        }

        auto range2 = scvfs(fvGeometry);
        NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        for (auto&& scvf : scvfs(fvGeometry))
        {
            std::cout << "-- scvf " << scvf.index() << " ip at: " << scvf.ipGlobal();
            if (scvf.boundary()) std::cout << " (on boundary).";
            std::cout << std::endl;
        }
    }
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {

    std::cout << e << std::endl;
    return 1;
}
