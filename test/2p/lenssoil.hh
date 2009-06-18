// $Id:$
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUNE_LENSSOIL_HH
#define DUNE_LENSSOIL_HH

/**
 * @file
 * @brief  Class for defining an instance of a Matrix2p soil
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class LensSoil: public Matrix2p<Grid,Scalar>
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::ctype DT;
    enum {dim=Grid::dimension, numEq=1};

    // define PERMEABILITY tensor
    virtual const FieldMatrix<DT,dim,dim> &K (const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
    {
        if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
            && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
            return Kin_;
        else
            return Kout_;
    }
    virtual double porosity(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
    {
        return 0.4;
    }

    virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T) const
    {
        if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
            && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
            return 0.18;
        else
            return 0.05;
    }

    virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T) const
    {
        return 0.0;
    }

    virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
    {
        return Matrix2p<Grid,Scalar>::van_genuchten;
    }

    virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T) const
    {
        // example for van Genuchten parameters
        std::vector<double> param(5);

        if ((x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0])
            && (x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1]))
        {
            param[0] = 1-1/4.7;
            param[1] = 4.7;
            param[2] = 0.5;
            param[3] = 1/3.;
            param[4] = 0.0037;
        }
        else
        {
            param[0] = 1-1/7.3;
            param[1] = 7.3;
            param[2] = 1/2.;
            param[3] = 1/3.;
            param[4] = 0.00045;

        }

        return param;
    }

    LensSoil() : Matrix2p<Grid,Scalar>()
    {
        Kin_ = Kout_ = 0;
        for(int i = 0; i < dim; i++)
        {
            Kin_[i][i] = 1e-13;
            Kout_[i][i] = 5e-10;
        }
    }

    ~LensSoil()
    {}

    //! Set the bounding box of the fine-sand lens
    void setLensCoords(const FieldVector<DT,dim>& innerLowerLeft,
                       const FieldVector<DT,dim>& innerUpperRight)
    {
        innerLowerLeft_ = innerLowerLeft;
        innerUpperRight_ = innerUpperRight;
    }

private:
    FieldMatrix<DT,dim,dim> Kin_;
    FieldMatrix<DT,dim,dim> Kout_;
    FieldVector<DT,dim> innerLowerLeft_;
    FieldVector<DT,dim> innerUpperRight_;
};

} // end namespace
#endif

