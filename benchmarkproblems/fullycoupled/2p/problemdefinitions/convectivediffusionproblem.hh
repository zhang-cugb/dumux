#ifndef DUNE_CONVECTIVEDIFFUSIONPROBLEM_HH
#define DUNE_CONVECTIVEDIFFUSIONPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>
#include<dumux/twophase/twophaseproblem_deprecated.hh>

/**
 * @file
 * @brief  Class defining a Convection-Diffusion problem
 * @author Markus Wolff
 */

namespace Dune {

template<class G, class RT> class ConvectiveDiffusionProblem :
	public TwoPhaseProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

public:
	enum {wPhaseIdx = 0, nwPhaseIdx = 1,pWIdx = 0, sNIdx = 1};
	enum {swrIdx = 0, snrIdx = 1};
	enum {VanGenuchten = 1, alphaIdx = 2, nIdx = 3};
	enum {BrooksCorey = 0,lambdaIdx = 2, p0Idx = 3};
	bool analytical_;

	virtual const FieldMatrix<DT,n,n>& K(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) {
		return K_;
	}

	virtual FieldVector<RT,m> q(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const {
		FieldVector<RT,m> values(0);

		return values;
	}

	virtual FieldVector<BoundaryConditions::Flags, m> bctype(
			const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const {
		FieldVector<BoundaryConditions::Flags, m> values(
				BoundaryConditions::neumann);

		if (x[0] < LowerLeft_[0] + eps_) {
			values[wPhaseIdx] = BoundaryConditions::dirichlet;
			values[nwPhaseIdx] = BoundaryConditions::dirichlet;
		}
//		if (x[0] > UpperRight_[0] - eps_) {
//			values[wPhaseIdx] = BoundaryConditions::dirichlet;
//			values[nwPhaseIdx] = BoundaryConditions::dirichlet;
//		}

		return values;
	}

	virtual FieldVector<RT,m> g(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const {
		FieldVector<RT,m> values(0);

		if (x[0] < LowerLeft_[0] + eps_) {
			values[pWIdx] = pwleftbc_;
			values[sNIdx] = 0.1;//Snr_;
		}
//		if (x[0] > UpperRight_[0] - eps_) {
//			values[pWIdx] = pwrightbc_;
//			values[sNIdx] = 1-Sinit_;
//		}

		return values;
	}

	virtual FieldVector<RT,m> J(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const {
		FieldVector<RT,m> values(0);

//		if (x[0] < LowerLeft_[0] + eps_) {
//			values[wPhaseIdx] = -0.0003;
//			values[nwPhaseIdx] = 0;
//		}

		if (x[0] > UpperRight_[0] - eps_) {
//			values[wPhaseIdx] = 0;
			values[nwPhaseIdx] = 0.0003;
		}

		return values;
	}

	virtual FieldVector<RT,m> initial(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const {

		FieldVector<RT,m> values(0);

		values[pWIdx] = pwleftbc_;
		values[sNIdx] = 1 - Sinit_ ;

		return values;
	}

	double porosity(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const {
		return Porosity_;
	}

	virtual FieldVector<RT,n> gravity() const {
		return gravity_;
	}

	virtual FieldVector<RT,4> materialLawParameters(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const {
		FieldVector<RT,4> values;

		if (chooselaw_) {
			values[swrIdx] = Swr_;
			values[snrIdx] = Snr_;
			values[alphaIdx] = Alpha_;
			values[nIdx] = N_;
		} else {
			values[swrIdx] = Swr_;
			values[snrIdx] = Snr_;
			values[lambdaIdx] = Lambda_;
			values[p0Idx] = p0_;
		}

		return values;
	}

	void assembleAdvPart(const G& g, const Entity& e, FieldVector<RT,m> flux);

	ConvectiveDiffusionProblem(TwoPhaseRelations& law = *(new LinearLaw), const FieldVector<DT,n> LowerLeft = 0,
			const FieldVector<DT,n> UpperRight = 0, int chooselaw = BrooksCorey,
			bool calcpeclet = false, bool analyticalsol = false, RT K = 1e-7,
			RT Swr = 0, RT Snr = 0, RT Si = 0, RT Porosity = 0.2,
			RT Lambda = 2.0, RT p0 = 10, RT Alpha = 0.1,//1.74e-4,
			RT N = 3.1257, RT pwleftbc=2e5, RT pwrightbc=2e5) :
		TwoPhaseProblem<G, RT>(law, analyticalsol),
				LowerLeft_(LowerLeft), UpperRight_(UpperRight), chooselaw_(chooselaw), K_(K),
				Swr_(Swr), Snr_(Snr), Sinit_(Si), Porosity_(Porosity),
				Lambda_(Lambda), p0_(p0), Alpha_(Alpha), N_(N),
				pwleftbc_(pwleftbc), pwrightbc_(pwrightbc),
				eps_(1e-8*UpperRight[0]), densityW_(law.wettingPhase.density()),
				densityN_(law.nonwettingPhase.density()) {
		switch (n) {
		case 1: //1D
			width_ = UpperRight[0] - LowerLeft[0];
			height_=1;
			gravity_[0]= 0;
			break;
		case 2: //2D
			K_[0][0]=K_[1][1]=K;
			K_[1][0]=K_[0][1]=0;

			height_ = UpperRight[1] - LowerLeft[1];
			width_ = UpperRight[0] - LowerLeft[0];

			gravity_[0] = 0;
			gravity_[1] = 0;
			break;
		default:
			DUNE_THROW(NotImplemented, "Dimension");
			break;
		}
	}

private:
	FieldMatrix<DT,n,n> K_;

	FieldVector<DT,n> LowerLeft_;
	FieldVector<DT,n> UpperRight_;

	DT width_, height_;
	DT eps_;
	RT densityW_, densityN_;
	FieldVector<DT,n> gravity_;
	RT Swr_, Snr_;
	RT Sinit_;
	RT Porosity_;
	RT Lambda_;
	RT p0_;
	RT Alpha_;
	RT N_;
	RT pwleftbc_, pwrightbc_;
	int chooselaw_;

};
}
#endif
