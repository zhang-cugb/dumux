// $Id$

#ifndef DUNE_TWOPHASEHEATPROBLEM_HH
#define DUNE_TWOPHASEHEATPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/property_baseclasses.hh>


/**
 * @file
 * @brief  Base class for defining an instance of the non-isothermal two-phase problem
 * @author Bernd Flemisch, Melanie Darcis
 */

namespace Dune {

template<class G, class RT> class TwoPhaseHeatProblem {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=3};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

public:
	//! evaluate source term
	/*! evaluate source term at given location
	 @param[in]  x    position in global coordinates
	 @param[in]  e    entity of codim 0
	 @param[in]  xi   position in reference element of e
	 \return     value of source term
	 */
	virtual FieldVector<RT,m> q(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const = 0;

	//! return type of boundary condition at the given global coordinate
	/*! return type of boundary condition at the given global coordinate
	 @param[in]  x    position in global coordinates
	 \return     boundary condition type given by enum in this class
	 */
	//	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,n>& x, const Entity& e,
	//			const IntersectionIterator& intersectionIt,
	//			const FieldVector<DT,n>& xi) const = 0;

	virtual FieldVector<BoundaryConditions::Flags, m>bctype(
			const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const = 0;

	//! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
		/*! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
		 @param[in]  x    position in global coordinates
		 \return     index of the primary variable
		 */

	virtual void dirichletIndex(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi, FieldVector<int,m>& dirichletIdx) const
	{
		for (int i = 0; i < m; i++)
			dirichletIdx[i]=i;
		return;
	}

	//! evaluate Dirichlet boundary condition at given position
	/*! evaluate Dirichlet boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,m> g(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const = 0;

	//! evaluate Neumann boundary condition at given position
	/*! evaluate Neumann boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,m> J(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const = 0;

	//! evaluate initial condition at given position
	/*! evaluate initial boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,m> initial(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const = 0;

	virtual FieldVector<RT,n> gravity() const = 0;

	//! properties of the wetting (liquid) phase
	/*! properties of the wetting (liquid) phase
	  \return	wetting phase
	 */
	virtual Fluid& wettingPhase () const
	{
		return wettingPhase_;
	}

	//! properties of the nonwetting (liquid) phase
	/*! properties of the nonwetting (liquid) phase
	  \return	nonwetting phase
	 */
	virtual Fluid& nonwettingPhase () const
	{
		return nonwettingPhase_;
	}

	//! properties of the soil
	/*! properties of the soil
	  \return	soil
	 */
    virtual Matrix2p<G, RT>& soil () const
    {
    	return soil_;
    }

	//! object for definition of material law
	/*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
	  \return	material law
	 */
    virtual TwoPhaseRelations<G, RT>& materialLaw () const
	{
		return materialLaw_;
	}

	TwoPhaseHeatProblem(Fluid& liq1, Fluid& liq2, Matrix2p<G, RT>& soil,
			TwoPhaseRelations<G,RT>& materialLaw = *(new TwoPhaseRelations<G,RT>))
	: wettingPhase_(liq1), nonwettingPhase_(liq2), soil_(soil),
	  materialLaw_(materialLaw)
	  { 	}

	//! always define virtual destructor in abstract base class
	virtual ~TwoPhaseHeatProblem() {
	}

protected:
	Fluid& wettingPhase_;
	Fluid& nonwettingPhase_;
    Matrix2p<G, RT>& soil_;
	TwoPhaseRelations<G, RT>& materialLaw_;
};

}
#endif
