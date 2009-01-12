// $Id$
#ifndef DUNE_VARIABLECLASS_HH
#define DUNE_VARIABLECLASS_HH

#include <dune/istl/bvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune
{

template<class Grid, class Scalar> class VariableClass
{
    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    typedef Dune::BlockVector< Dune::FieldVector<Scalar,dim> > SlopeType;

typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::Traits::LevelIndexSet IndexSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout>
    ElementMapper;

public:
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarVectorType;
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;

    Grid& grid;
    int transLevel;
    int diffLevel;
    ElementMapper diffMapper;
    ElementMapper transMapper;
    int diffSize;
    int transSize;

    ScalarVectorType saturation;
    ScalarVectorType pressure;
    VelType velocity;
    SlopeType slope;

    VariableClass(Grid& grid, Scalar& initialSat = *(new Scalar(0)), Scalar& initalPress = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)), int transLev = -1, int diffLev = -1)
    : grid(grid),
    transLevel((transLev >= 0) ? transLev : grid.maxLevel()), diffLevel((diffLev >= 0) ? diffLev : grid.maxLevel()),
    diffMapper(grid, (diffLev >= 0) ? grid.levelIndexSet(diffLev) : grid.levelIndexSet(grid.maxLevel())), transMapper(grid, (transLev >= 0) ? grid.levelIndexSet(transLev) : grid.levelIndexSet(grid.maxLevel())),
    diffSize(diffMapper.size()),transSize(transMapper.size())
    {
        initSat(initialSat, transSize);
        initPress(initalPress, diffSize);
        initVel(initialVel, transSize);
        initSlopes(0.0, transSize);
    }

    void initSat(Scalar& initialSat, int size)
    {
        saturation.resize(size);
        saturation=initialSat;
        return;
    }

    void initSlopes(Scalar initialSlope, int size)
    {
        slope.resize(size);
        slope=initialSlope;
        return;
    }

    void initPress(Scalar& initialPress, int size)
    {
        pressure.resize(size);
        pressure=initialPress;
        return;
    }
    void initVel(Dune::FieldVector<Scalar, dim>& initialVel, int size)
    {
        velocity.resize(size);
        velocity=initialVel;
        return;
    }

    ScalarVectorType& sat() const
    {
        return saturation;
    }
    ScalarVectorType& press() const
    {
        return pressure;
    }
    VelType& vel() const
    {
        return velocity;
    }

    const Dune::FieldVector<Scalar,1>& sat(const GlobalPosition globalPos,
            const Element& element, const LocalPosition localPos) const
    {
        return saturation[transMapper.map(element)];;
    }

    const Dune::FieldVector<Scalar,1>& press(const GlobalPosition globalPos,
            const Element& element, const LocalPosition localPos) const
    {
        return pressure[diffMapper.map(element)];
    }

    const Dune::FieldVector<Scalar,dim>& vTotal(const Element& element,
            const int numberInSelf) const
    {
        int elemId = transMapper.map(element);

        return (velocity[elemId][numberInSelf]);
    }

    /*! @brief prints the saturation to a VTK file
     *
     *  The file name is "<name>-<k>.vtu" where k is an integer number.
     *  @param name specifies the name of the VTK file
     *  @param k specifies a number
     */

    void vtkout (const char* name, int k)
    {
        vtkoutMultiLevel(name,k,diffLevel,transLevel);
    }

    void vtkoutMultiLevel(const char* name, int k, int pressureLevel = 0, int satLevel = 0) const
    {
        if (pressureLevel == satLevel)
        {
            VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
            char fname[128];
            sprintf(fname, "%s-%05d", name, k);
            vtkwriter.addCellData(saturation, "saturation");
            vtkwriter.addCellData(pressure, "total pressure p~");
            vtkwriter.write(fname, VTKOptions::ascii);
        }
        else
        {
            Dune::VTKWriter<typename Grid::LevelGridView>
            vtkwriterpressure(grid.levelView(pressureLevel));
            char fname[128];
            sprintf(fname, "%s-press%05d", name, k);
            vtkwriterpressure.addCellData(pressure, "total pressure p~");
            vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

            Dune::VTKWriter<typename Grid::LevelGridView>
            vtkwritersaturation(grid.levelView(satLevel));
            sprintf(fname, "%s-%05d", name, k);
            vtkwritersaturation.addCellData(saturation, "saturation");
            vtkwritersaturation.write(fname, VTKOptions::ascii);
        }
        return;
    }
    //    void vtkoutpressure(const char* name, int k) const
    //    {
    //        LeafVTKWriter<typename Grid::LeafGridView> vtkwriter(grid);
    //        char fname[128];
    //        sprintf(fname, "%s-press%05d", name, k);
    //        vtkwriter.addCellData(pressure, "total pressure p~");
    //        vtkwriter.write(fname, VTKOptions::ascii);
    //    }
    //    void vtkoutsamelevel(const char* name, int k, BlockVector<FieldVector<Scalar, 2> > uEx) const
    //    {
    //        ScalarVectorType saturationEx(transsize);
    //        ScalarVectorType error(transsize);
    //
    //        for (int i=0; i<transsize; i++)
    //        {
    //            saturationEx[i] = uEx[i][0];
    //            error[i]=uEx[i][1];
    //        }
    //        VTKWriter<typename Grid::LevelGridView> vtkwriter(grid);
    //        char fname[128];
    //        sprintf(fname, "%s-%05d", name, k);
    //        vtkwriter.addCellData(saturation, "saturation");
    //        vtkwriter.addCellData(pressure, "total pressure p~");
    //        vtkwriter.addCellData(saturationEx, "saturation (exact solution)");
    //        vtkwriter.addCellData(error, "error");
    //        vtkwriter.write(fname, VTKOptions::ascii);
    //    }
};
}
#endif
