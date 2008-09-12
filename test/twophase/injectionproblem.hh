#ifndef DUNE_INJECTIONPROBLEM_HH
#define DUNE_INJECTIONPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/multicomponentrelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/2p2c/2p2cproblem.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhaseTwoComponent problem
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
  //! base class that defines the parameters of a diffusion equation
  /*! An interface for defining parameters for the stationary diffusion equation
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
   * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
   * on \f$\Gamma_2\f$. Here,
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
   * and \f$\lambda\f$ the total mobility, possibly depending on the
   * saturation.
   *
   *	Template parameters are:
   *
   *	- Grid  a DUNE grid type
   *	- RT    type used for return values
   */
  template<class G, class RT>
  class InjectionProblem : public TwoPTwoCProblem<G, RT>
  {
	  enum {dim=G::dimension, m=2};
	  typedef typename G::ctype DT;
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	  enum {pWIdx = 0, switchIdx = 1}; // phase index
	  enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // phase state


/////////////////////////////
// TYPE of the boundaries
/////////////////////////////
	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<DT,dim>& xi) const
	{
		FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

		if (x[0] < outerLowerLeft_[0] + eps_)
			values = BoundaryConditions::dirichlet;
		if (x[1] < eps_)
			values = BoundaryConditions::dirichlet;

		return values;
	}

/////////////////////////////
// DIRICHLET boundaries
/////////////////////////////
	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		values[pWIdx] = 1e5 - densityW_*gravity_[1]*(depthBOR_ - x[1]);
		values[switchIdx] = 0.05;

//		if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//		 && x[0] >= innerLowerLeft_[0])
//			values[switchIdx] = 0.2;
//		else
//			values[switchIdx] = 1e-6;

		return values;
	}

/////////////////////////////
// NEUMANN boundaries
/////////////////////////////
	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		//RT lambda = (x[1])/height_;

		if (x[1] < 15.0 && x[1] > 5.0)
			values[switchIdx] = -1e-2;

		return values;
	}

/////////////////////////////
// INITIAL values
/////////////////////////////
	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
				const FieldVector<DT,dim>& xi) const
	{

		FieldVector<RT,m> values;

		values[pWIdx] = 1e5 - densityW_*gravity_[1]*(depthBOR_ - x[1]);
		values[switchIdx] = 0.01;

//			if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//			 && x[0] >= innerLowerLeft_[0])
//				values[switchIdx] = 0.2;
//			else
//				values[switchIdx] = 1e-6;

		return values;
	}


	int initialPhaseState (const FieldVector<DT,dim>& x, const Entity& e,
				  const FieldVector<DT,dim>& xi) const
	{

		enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
		int state;

//			if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//			      && x[0] >= innerLowerLeft_[0])
//				state = 2;
//			else

		state = bothPhases;

		return state;
	}

/////////////////////////////
// sources and sinks
/////////////////////////////
	virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e,
					const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

//////////////////////////////

	virtual FieldVector<RT,dim> gravity () const
	{
		return gravity_;
	}

	double depthBOR () const
	{
		return depthBOR_;
	}

	InjectionProblem(Liquid_GL& liq, Gas_GL& gas, Matrix2p<G, RT>& soil,
			const FieldVector<DT,dim> outerLowerLeft = 0., const FieldVector<DT,dim> outerUpperRight = 0.,
			const FieldVector<DT,dim> innerLowerLeft = 0., const FieldVector<DT,dim> innerUpperRight = 0.,
			const RT depthBOR = 0., TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>),
			MultiComp& multicomp = *(new CWaterAir))
	: TwoPTwoCProblem<G,RT>(liq, gas, soil, multicomp, law),
	  outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
	  depthBOR_(depthBOR), eps_(1e-8*outerUpperRight[0])
	  {
		height_ = outerUpperRight[1] - outerLowerLeft[1];
		width_ = outerUpperRight[0] - outerLowerLeft[0];

		gravity_[0] = 0;
		gravity_[1] = -9.81;
	  }

  private:
	  FieldVector<DT,dim> outerLowerLeft_, outerUpperRight_;
	  FieldVector<DT,dim> innerLowerLeft_, innerUpperRight_;
	  DT width_, height_;
	  DT depthBOR_, eps_;
	  RT densityW_, densityN_;
	  FieldVector<DT,dim> gravity_;
  };

 ////////////////////////////////////////--SOIL--//////////////////////////////////////////////
  template<class G, class RT>
  class InjectionSoil: public Matrix2p<G,RT>
  {
  public:
  	typedef typename G::Traits::template Codim<0>::Entity Entity;
  	typedef typename G::ctype DT;
  	enum {dim=G::dimension, m=1};

  	virtual FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
  	{
  		return K_;
  	}
  	virtual double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return 0.3;
  	}

  	virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		return 0.2;
  	}

  	virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		return 0.05;
  	}

  	/* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
  			 * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
  	virtual double heatCap(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return 	790 /* spec. heat cap. of granite */
  						* 2700 /* density of granite */
  						* porosity(x, e, xi);
  	}

  	virtual double heatCond(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double sat) const
  	{
  		static const double lWater = 0.6;
  		static const double lGranite = 2.8;
  		double poro = porosity(x, e, xi);
  		double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
  		double ldry = pow(lGranite, (1-poro));
  		return ldry + sqrt(sat) * (ldry - lsat);
  	}

  	virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
  	{
  		// example for Brooks-Corey parameters
  		std::vector<double> param(2);
  		param[0] = 2.; // lambda
  		param[1] = 1e4; // entry-pressures

  		return param;
  	}

  	virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
  	{
  		return Matrix2p<G,RT>::brooks_corey;
  	}

  	InjectionSoil():Matrix2p<G,RT>()
  	{
  		for(int i = 0; i < dim; i++)
  			K_[i][i] = 1e-12;
  	}

  	~InjectionSoil()
  	{}

  private:
  	FieldMatrix<DT,dim,dim> K_;

  };


} //end namespace
#endif
