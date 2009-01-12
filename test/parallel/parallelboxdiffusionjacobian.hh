#ifndef DUNE_PARALLELBOXDIFFUSIONJACOBIAN_HH
#define DUNE_PARALLELBOXDIFFUSIONJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/operators/boxjacobian.hh>
#include "diffusionparameters.hh"

namespace Dune
{
  //! A class for computing local jacobian matrices
  /*! A class for computing local jacobian matrix for the
    diffusion equation

        div j = q; j = -K grad u; in Omega

        u = g on Gamma1; j*n = J on Gamma2.

    Uses the box method.

    Template parameters are:

    - G     a DUNE grid type
    - RT    type used for return values
  */
  template<class G, class RT, class BoxFunction = LeafP1FunctionExtended<G, RT, 1> >
  class ParallelBoxDiffusionJacobian
    : public BoxJacobian<ParallelBoxDiffusionJacobian<G,RT,BoxFunction>,G,RT,1,BoxFunction>
  {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef ParallelBoxDiffusionJacobian<G,RT,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::VBlockType VBlockType;
     typedef FVElementGeometry<G> FVElementGeometry;

  public:
    enum {n=G::dimension};

    //! Constructor
    ParallelBoxDiffusionJacobian (DiffusionParameters<G,RT>& params,
                  bool levelBoundaryAsDirichlet_, const G& grid,
                  BoxFunction& sol,
                  bool procBoundaryAsDirichlet_=true)
    : BoxJacobian<ThisType,G,RT,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
      problem(params)
    {
      this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
        VBlockType result(0);
        return result;
    }

    VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
        return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
    }

    VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
        FieldVector<RT, n> gradP(0);
        for (int k = 0; k < this->fvGeom.numVertices; k++) {
            FieldVector<RT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k];
            gradP += grad;
        }

        FieldVector<RT,n> KGradP(0);
        elData.K.umv(gradP, KGradP);

        VBlockType flux = KGradP*this->fvGeom.subContVolFace[face].normal;

        return flux;
    }

    void computeElementData (const Entity& e)
    {
        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        elData.K = problem.K(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
    };

    void updateStaticData (const Entity& e, const VBlockType* sol)
    {
        return;
    }

    void updateVariableData (const Entity& e, const VBlockType* sol, int i, bool old = false)
    {
        return;
    }

    void updateVariableData (const Entity& e, const VBlockType* sol, bool old = false)
    {
        return;
    }

    struct ElementData {
        FieldMatrix<DT,n,n> K;
       };

       ElementData elData;
    DiffusionParameters<G,RT>& problem;
  };
}
#endif
