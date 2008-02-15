#ifndef DUNE_FVELEMENTGEOMETRY_HH
#define DUNE_FVELEMENTGEOMETRY_HH

namespace Dune
{

template<class G>
class FVElementGeometry  
{
	enum{dim = G::dimension};
	enum{maxNC = (dim < 3 ? 4 : 8)};
	typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;

	struct SDValues 
	{
		FieldVector<DT, maxNC> shape;                    /* values of shape functions at ip              */
		FieldVector<FieldVector<DT, dim>, maxNC> grad;              /* derivatives of shape functions at ip */
		FieldMatrix<DT, dim, dim> J;                         /* jacobian at ip                               */
		FieldMatrix<DT, dim, dim> Jinv;                  /* inverse of jacobian at ip                    */
		DT detJ;                                    /* det of jacobian at ip                                */
	};

	struct SubControlVolume /* FV intersected with element                  */
	{
		int co;                                                 /* # of corner                                                  */
		FieldVector<DT, dim> center;                   /* node position                                                */
	    DT volume;                                  /* volume (area) of scv                                 */
	};                                     

	struct SubControlVolumeFace
	{
		int i,j;                                                /* scvf seperates corner i and j of elem*/
		FieldVector<DT, dim> ip_local;                 /* integration point in local coords    */
		FieldVector<DT, dim> ip_global;                    /* integration point in global coords       */
		FieldVector<DT, dim> normal;                   /* normal on face at ip pointing to CV j*/
		SDValues sdv;                                  /* shape fcts, deriv. etc. at scv-faces */
	};

	typedef struct {
		int co;                                                 /* corresponding corner                                 */
		int side;                                               /* boundary side of element                             */
		FieldVector<DT, dim> ip_local;                 /* integration point in local coords    */
		FieldVector<DT, dim-1> param;            /* local side coordinates                       */
		FieldVector<DT, dim> normal;                   /* normal on face at ip pointing to CV j*/
		DT area;                                    /* area of boundary face                                */
		SDValues sdv;                                  /* shape fcts, deriv. etc. at b-faces   */
	} BoundaryFace; 

};

}


#endif



