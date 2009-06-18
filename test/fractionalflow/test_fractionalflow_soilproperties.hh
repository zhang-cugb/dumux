// $Id$
#ifndef TEST_FRACTIONALFLOW_SOILPROPERTIES_HH
#define TEST_FRACTIONALFLOW_SOILPROPERTIES_HH

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class FractionalFlowTestSoil: public Matrix2p<Grid, Scalar>
{
public:
    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Matrix2p<Grid, Scalar>::modelFlag modelFlag;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    virtual const FieldMatrix &K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return constPermeability_;
    }
    virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return 0.2;
    }

    virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.0;
    }

    virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0.0;
    }


    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T = 283.15) const
    {
        return brooksCoreyParameters_;
    }

    virtual modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return Matrix2p<Grid, Scalar>::brooks_corey;
    }

    FractionalFlowTestSoil(Scalar entryPressure = 0)
    : Matrix2p<Grid,Scalar>(),constPermeability_(0), entryPressure_(entryPressure),
    brooksCoreyParameters_(2)
    {
        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-7;
        }

        brooksCoreyParameters_[0] = 2.0; // lambda
        brooksCoreyParameters_[1] = entryPressure_;
    }

private:
    FieldMatrix constPermeability_;
    Scalar entryPressure_;
    std::vector<double> brooksCoreyParameters_;

public:

    //    RandomPermeability<Grid> permeability;
};

} // end namespace
#endif /*MATRIXPROPERTIES*/
