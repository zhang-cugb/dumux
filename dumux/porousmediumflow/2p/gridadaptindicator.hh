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
 * \ingroup TwoPModel
 * \brief Class defining a standard, saturation dependent indicator for grid adaptation
 */

#ifndef DUMUX_TWOP_ADAPTION_INDICATOR_HH
#define DUMUX_TWOP_ADAPTION_INDICATOR_HH

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux
{

/*!\ingroup TwoPModel
 * \brief  Class defining a standard, saturation dependent indicator for grid adaptation
 */
template<class TypeTag>
class TwoPGridAdaptIndicator
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    enum { saturationIdx = Indices::saturationIdx };

public:
    /*! \brief The Constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     *
     *  Note: refineBound_, coarsenBound_ & maxSaturationDelta_ are chosen
     *        in a way such that the indicator returns false for all elements
     *        before having been calculated.
     */
    TwoPGridAdaptIndicator(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : fvGridGeometry_(fvGridGeometry)
    , refineBound_(std::numeric_limits<Scalar>::max())
    , coarsenBound_(std::numeric_limits<Scalar>::lowest())
    , maxSaturationDelta_(fvGridGeometry_->gridView().size(0), 0.0)
    , minLevel_(getParamFromGroup<std::size_t>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Adaptive.MinLevel", 0))
    , maxLevel_(getParamFromGroup<std::size_t>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Adaptive.MaxLevel", 0))
    {}

    /*!
     * \brief Function to set the minimum allowed level.
     *
     * \param minLevel The minimum level
     */
    void setMinLevel(std::size_t minLevel)
    {
        minLevel_ = minLevel;
    }

    /*!
     * \brief Function to set the maximum allowed level.
     *
     *\param maxLevel The maximum level
     */
    void setMaxLevel(std::size_t maxLevel)
    {
        maxLevel_ = maxLevel;
    }

    /*!
     * \brief Function to set the minumum/maximum allowed levels.
     *
     * \param minLevel The minimum level
     * \param maxLevel The maximum level
     */
    void setLevels(std::size_t minLevel, std::size_t maxLevel)
    {
        minLevel_ = minLevel;
        maxLevel_ = maxLevel;
    }

    /*!
     * \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     * \param sol The solution vector
     * \param refineTol The refinement tolerance
     * \param coarsenTol The coarsening tolerance
     *
     *  This standard two-phase indicator is based on the saturation gradient.
     */
    void calculate(const SolutionVector& sol,
                   Scalar refineTol = 0.05,
                   Scalar coarsenTol = 0.001)
    {
        //! reset the indicator to a state that returns false for all elements
        refineBound_ = std::numeric_limits<Scalar>::max();
        coarsenBound_ = std::numeric_limits<Scalar>::lowest();
        maxSaturationDelta_.assign(fvGridGeometry_->gridView().size(0), 0.0);

        //! maxLevel_ must be higher than minLevel_ to allow for refinement
        if (minLevel_ >= maxLevel_)
            return;

        //! check for inadmissible tolerance combination
        if (coarsenTol > refineTol)
            DUNE_THROW(Dune::InvalidStateException, "Refine tolerance must be higher than coarsen tolerance");

        //! variables to hold the max/mon saturation values on the leaf
        Scalar globalMax = std::numeric_limits<Scalar>::lowest();
        Scalar globalMin = std::numeric_limits<Scalar>::max();

        //! Calculate minimum and maximum saturation
        for (const auto& element : elements(fvGridGeometry_->gridView()))
        {
            //! index of the current leaf-element
            const auto globalIdxI = fvGridGeometry_->elementMapper().index(element);

            //! obtain the saturation at the center of the element
            const auto geometry = element.geometry();
            const ElementSolution elemSol(element, sol, *fvGridGeometry_);
            const Scalar satI = evalSolution(element, geometry, *fvGridGeometry_, elemSol, geometry.center())[saturationIdx];

            //! maybe update the global minimum/maximum
            using std::min;
            using std::max;
            globalMin = min(satI, globalMin);
            globalMax = max(satI, globalMax);

            //! calculate maximum delta in saturation for this cell
            for (const auto& intersection : intersections(fvGridGeometry_->gridView(), element))
            {
                //! Only consider internal intersections
                if (intersection.neighbor())
                {
                    //! Access neighbor
                    const auto outside = intersection.outside();
                    const auto globalIdxJ = fvGridGeometry_->elementMapper().index(outside);

                    //! Visit intersection only once
                    if (element.level() > outside.level() || (element.level() == outside.level() && globalIdxI < globalIdxJ))
                    {
                        //! obtain saturation in the neighbor
                        const auto outsideGeometry = outside.geometry();
                        const ElementSolution elemSolJ(outside, sol, *fvGridGeometry_);
                        const Scalar satJ = evalSolution(outside, outsideGeometry, *fvGridGeometry_, elemSolJ, outsideGeometry.center())[saturationIdx];

                        using std::abs;
                        Scalar localdelta = abs(satI - satJ);
                        maxSaturationDelta_[globalIdxI] = max(maxSaturationDelta_[globalIdxI], localdelta);
                        maxSaturationDelta_[globalIdxJ] = max(maxSaturationDelta_[globalIdxJ], localdelta);
                    }
                }
            }
        }

        //! compute the maximum delta in saturation
        const auto globalDelta = globalMax - globalMin;

        //! compute the refinement/coarsening bounds
        refineBound_ = refineTol*globalDelta;
        coarsenBound_ = coarsenTol*globalDelta;

// TODO: fix adaptive simulations in parallel
//#if HAVE_MPI
//    // communicate updated values
//    using DataHandle = VectorExchange<ElementMapper, ScalarSolutionType>;
//    DataHandle dataHandle(problem_.elementMapper(), maxSaturationDelta_);
//    problem_.gridView().template communicate<DataHandle>(dataHandle,
//                                                         Dune::InteriorBorder_All_Interface,
//                                                         Dune::ForwardCommunication);
//
//    using std::max;
//    refineBound_ = problem_.gridView().comm().max(refineBound_);
//    coarsenBound_ = problem_.gridView().comm().max(coarsenBound_);
//
//#endif

        //! check if neighbors have to be refined too
        for (const auto& element : elements(fvGridGeometry_->gridView(), Dune::Partitions::interior))
            if (this->operator()(element) > 0)
                checkNeighborsRefine_(element);
    }

    /*! \brief function call operator to return mark
     *
     *  \return  1 if an element should be refined
     *          -1 if an element should be coarsened
     *           0 otherwise
     *
     *  \param element A grid element
     */
    int operator() (const Element& element) const
    {
        if (element.hasFather()
            && maxSaturationDelta_[fvGridGeometry_->elementMapper().index(element)] < coarsenBound_)
        {
            return -1;
        }
        else if (element.level() < maxLevel_
                 && maxSaturationDelta_[fvGridGeometry_->elementMapper().index(element)] > refineBound_)
        {
            return 1;
        }
        else
            return 0;
    }

private:
    /*!
     * \brief Method ensuring the refinement ratio of 2:1
     *
     *  For any given element, a loop over the neighbors checks if the
     *  entities refinement would require that any of the neighbors has
     *  to be refined, too. This is done recursively over all levels of the grid.
     *
     * \param element Element of interest that is to be refined
     * \param level level of the refined element: it is at least 1
     * \return true if everything was successful
     */
    bool checkNeighborsRefine_(const Element &element, std::size_t level = 1)
    {
        for(const auto& intersection : intersections(fvGridGeometry_->gridView(), element))
        {
            if(!intersection.neighbor())
                continue;

            // obtain outside element
            const auto outside = intersection.outside();

            // only mark non-ghost elements
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if (outside.level() < maxLevel_ && outside.level() < element.level())
            {
                // ensure refinement for outside element
                maxSaturationDelta_[fvGridGeometry_->elementMapper().index(outside)] = std::numeric_limits<Scalar>::max();
                if(level < maxLevel_)
                    checkNeighborsRefine_(outside, ++level);
            }
        }

        return true;
    }

    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    Scalar refineBound_;
    Scalar coarsenBound_;
    std::vector< Scalar > maxSaturationDelta_;
    std::size_t minLevel_;
    std::size_t maxLevel_;
};

} // end namespace Dumux

#endif /* DUMUX_TWOP_ADAPTION_INDICATOR_HH */
